#include <fstream>
#include "RigidBodyDemo.h"

#include "matplotlibcpp.h"
namespace plt = matplotlibcpp;

void saveMatrixToGNUPlotFile(const std::string& filename, const MatrixXd& matrix) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }

    int numRows = matrix.rows();
    int numCols = matrix.cols();

    for (int i = 0; i < numRows; ++i) {
        for (int j = 0; j < numCols; ++j) {
            file << j << " " << matrix(i, j) << "\n";
        }
        file << "\n";  // Separate rows by a blank line
    }

    file.close();
}

inline double W_fischer(const MatrixXd& D, const VectorXd& qdot, const VectorXd& lambda_f,
                        const double& r, const double& mu, const double& lambda_n)
{
    double ret = 1;
    if (lambda_n * mu > 0)
    {
        ret = r*(sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)) - r*(mu*lambda_n - lambda_f.norm()))/
               ((D.transpose()*qdot).squaredNorm()+r*mu*lambda_n-sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)));
    }
    if (isnan(ret))
        ret = 0;

    return ret;
}

inline double fischer(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    return a + b - sqrt(pow(a,2) + pow(b,2)) ;
}

inline double fischer_partial_a(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    double ret;
    if (a == 0 && b == 0)
    {
        ret = 0;
    }
    else
    {
        ret = 1 - (a/sqrt(pow(a,2) + pow(b,2)));
    }
    return ret;
}

inline double fischer_partial_b(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    double ret;
    if (a == 0 && b == 0)
    {
        ret = 1;
    }
    else
    {
        ret = 1 - (b/sqrt(pow(a,2) + pow(b,2)));
    }
    return ret;
}

VectorXd RigidBodyDemo::computeResidual(const VectorXd &x, const Quaterniond& q , const VectorXd &u,
                                    const VectorXd &u_tilde, const VectorXd& lambda,  const MatrixXd &M_tilde,
                                    const vector<Contact*>& contacts, double dt) {
    int num_contacts = contacts.size();
    VectorXd r(6 + lambda.rows());
    VectorXd phi_normal(num_contacts);
    VectorXd phi_friction(2*num_contacts);
    MatrixXd J_normal(num_contacts,6);
    MatrixXd J_friction(2*num_contacts,6);
    VectorXd lambda_normal = lambda.head(num_contacts);
    VectorXd lambda_friction = lambda.tail(2*num_contacts);

    for (int cidx = 0; cidx < num_contacts; cidx++) {

        // preconditioning values
        double r_n = 1; double r_f = 1;

        // lambda values
        lambda_normal = lambda.head(num_contacts);
        lambda_friction = lambda.tail(2*num_contacts);

        // the collision
        auto &contact = contacts[cidx];

        // the contact basis
        MatrixXd n_t_b = contact->n_t_b;

        // the normal
        Vector3d n = n_t_b.col(0);
        MatrixXd t_b = n_t_b.block<3,2>(0,1);

        // object indices
        int objectA_idx = contact->objectA->m_idx;
        int objectB_idx = contact->objectB->m_idx;

        // fill the normal constraint vector and normal jacobian
        if (objectA_idx == 0 )
        {
            // normal constraints
            double C_n = n.dot( ( x + q.toRotationMatrix()*contact->contact_wrt_objectA ) -
                                ( m_objects[1]->x + m_objects[1]->q.toRotationMatrix()*contact->contact_wrt_objectB ) ) - 1e-6;
            if (C_n > 0)
            {
                contact->m_is_active = false;
            }
            phi_normal(cidx) = fischer(C_n,lambda_normal(cidx)) / dt;
            double pfpa = fischer_partial_a(C_n, r_n * lambda_normal(cidx));
            double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal(cidx));
            J_normal.block<1, 3>(cidx, 0) = pfpa*n;
            J_normal.block<1, 3>(cidx, 3) = -pfpa*n.transpose()*vec2skew(contact->contact_wrt_objectA);

            // friction constraints
            Matrix<double,6,2> D; D.block<3,2>(0,0) = t_b;
            D.block<3,1>(3,0) = -t_b.block<3,1>(0,0).transpose()*vec2skew(contact->contact_wrt_objectA);
            D.block<3,1>(3,1) = -t_b.block<3,1>(0,1).transpose()*vec2skew(contact->contact_wrt_objectA);
            J_friction.block<2, 6>(2 * cidx, 0) = D.transpose();
            double W_f = W_fischer(D,u, lambda_friction.block<2,1>(2*cidx,0), r_f, 0.9, lambda_normal(cidx));
            phi_friction.block<2,1>(2*cidx,0) = D.transpose() * u + W_f * lambda_friction.block<2,1>(2*cidx,0);
        }
        else if (objectA_idx == 1)
        {
            // normal constraints
            double C_n = n.dot(( m_objects[1]->x + m_objects[1]->q.toRotationMatrix()*contact->contact_wrt_objectA ) -
                               ( x + q.toRotationMatrix()*contact->contact_wrt_objectB ) ) - 1e-6;

            if (C_n > 0)
            {
                contact->m_is_active = false;
            }

            phi_normal(cidx) = fischer(C_n,lambda_normal(cidx)) / dt;
            double pfpa = fischer_partial_a(C_n, r_n * lambda_normal(cidx));
            double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal(cidx));
            J_normal.block<1, 3>(cidx, 0) = -pfpa*n;
            J_normal.block<1, 3>(cidx, 3) = pfpa*n.transpose()* vec2skew(contact->contact_wrt_objectB);

            // friction constraints
            Matrix<double,6,2> D; D.block<3,2>(0,0) = -t_b;
            D.block<3,1>(3,0) = t_b.block<3,1>(0,0).transpose()*vec2skew(contact->contact_wrt_objectB);
            D.block<3,1>(3,1) = t_b.block<3,1>(0,1).transpose()*vec2skew(contact->contact_wrt_objectB);
            J_friction.block<2, 6>( 2 * cidx, 0) = D.transpose();
            double W_f = W_fischer(D,u, lambda_friction.block<2,1>(2*cidx,0), r_f, 0.9, lambda_normal(cidx));
            phi_friction.block<2,1>(2*cidx,0) = D.transpose() * u + W_f * lambda_friction.block<2,1>(2*cidx,0);
        }

    }

    r.block<6,1>(0,0) = M_tilde*(u - u_tilde)/dt - J_normal.transpose()*lambda_normal - J_friction.transpose()*lambda_friction;
    r.block(6,0,num_contacts,1) = phi_normal;
    r.block(6+num_contacts,0,2*num_contacts,1) = phi_friction;

    return r;
}

void RigidBodyDemo::backtrackingLineSearch(double &t, double &err,
                                       const VectorXd &x, const VectorXd &u,
                                       const VectorXd &u_tilde, const VectorXd &delta_u,
                                       const VectorXd& lambda, const VectorXd& delta_lambda, const MatrixXd &M_tilde,
                                       const vector<Contact*>& contacts, double dt, double alpha, double beta, int max_iter) {
//    VectorXd grad_err = delta_lambda_n;
//    double err_k;
//
//    for (int k = 0; k < max_iter; ++k) {
//        VectorXd lambda_n_k = lambda_n + t * delta_lambda_n;
//        VectorXd u_k = u + t * delta_u;
//        VectorXd x_k = x + u * dt;
//        VectorXd r = computeResidual(x_k, u_k, u_tilde, lambda_n_k, M_tilde, colInfo, dt);
//        err_k = 0.5 * r.dot(r);
//
//        // Perform Armijo condition to check if we have sufficient decrease
//        if (err_k <= err + beta * t * grad_err.dot(delta_lambda_n)) {
//            break;
//        }
//
//        t = alpha * t;
//    }
//
//    err = err_k;
}

void RigidBodyDemo::initialize()
{
    // add balls to the scene
    for (int i = 0 ; i < 2; i++)
    {
        if (i == 0)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/aldo/ImplicitNonlinearComplementarity/resources/PegInHole/peg.obj"));
            m_objects.back()->set_local_pos(Vector3d(0.1,0.2,1));
            auto v = m_objects.back()->m_vertices;
            Vector3d sum(0,0,0);
            for (int k = 0 ; k < v.rows();k++)
            {
                sum += v.row(k);
            }
            sum = sum /  v.rows();
            cout << "centroid = " << sum.transpose() << endl;
            m_objects.back()->scaleMesh(0.025);
            m_objects.back()->setWireMode(true);
            m_objects.back()->getMesh(0)->m_material->setBlue();
            m_objects.back()->setWireMode(true);
        }
        else if (i == 1)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/aldo/ImplicitNonlinearComplementarity/resources/PegInHole/block.obj"));
            m_objects.back()->set_is_static(true);
            auto v = m_objects.back()->m_vertices;
            Vector3d sum(0,0,0);
            for (int k = 0 ; k < v.rows();k++)
            {
                sum += v.row(k);
            }
            sum = sum /  v.rows();
            cout << "centroid = " << sum.transpose() << endl;
            m_objects.back()->scaleMesh(0.025);
            m_objects.back()->getMesh(0)->m_material->setBrownCornsilk();
            m_objects.back()->setWireMode(true);
        }

        m_world->addChild(m_objects.back());
        m_objects.back()->compute_inertia_matrix();
    }

    // create the collision detector
    m_collisionDetector = new CollisionDetector(m_objects);
}

void RigidBodyDemo::step(double dt)
{
    // reset forces
    for(auto b : m_objects)
    {
        b->f = b->mass * Vector3d(0., 0, -0.098);
        b->tau = Vector3d(0.0,0,0);
        b->fc.setZero();
        b->tauc.setZero();
        b->m_contacts.clear();
    }

    // computes any collisions
    m_collisionDetector->computeCollisions();

    // computes the constraints
    if (m_collisionDetector->m_contacts.size() > 0)
    {
        computeConstraints(dt);
    }
    else
    {
        m_objects[0]->xdot = m_objects[0]->xdot_tilde;
        m_objects[0]->x = m_objects[0]->x_tilde;
        m_objects[0]->omega = m_objects[0]->omega_tilde;
        m_objects[0]->q = m_objects[0]->q_tilde;
    }

    // explicit step
    m_objects[0]->update_inertia_matrix();
    m_objects[0]->xdot_tilde = m_objects[0]->xdot + dt * (1./m_objects[0]->mass) * (m_objects[0]->f);
    m_objects[0]->x_tilde = m_objects[0]->x + dt * m_objects[0]->xdot_tilde;
    m_objects[0]->omega_tilde = m_objects[0]->omega + dt * m_objects[0]->Iinv * (m_objects[0]->tau + m_objects[0]->tauc);
    m_objects[0]->q_tilde = sumQuaternions(m_objects[0]->q,
                          multiplyQuaternionScalar(dt * 0.5 , angularVelocityToQuaternion(m_objects[0]->omega) * m_objects[0]->q));
    m_objects[0]->q.normalize();

    // clears the contacts
    m_collisionDetector->clear();

}

void RigidBodyDemo::updateGraphics()
{
    for (auto object : m_objects)
    {
        object->update_mesh_position();
    }
}

void RigidBodyDemo::computeConstraints(double dt)
{
    auto contacts = m_collisionDetector->m_contacts; // get contacts
    int num_objects = m_objects.size(); // get the number of objects
    int num_contacts = contacts.size(); // get the number of contatcs
    cout << "Num Contacts = " << num_contacts << endl;

    VectorXd lambda(num_contacts*3); lambda.setZero(); // the generalized constraint vector
    VectorXd h(num_contacts*3); h.setZero();
    MatrixXd J(num_contacts*3,6); J.setZero();
    MatrixXd C(num_contacts*3,num_contacts*3); C.setZero();

    // fill the unconstrained positions
    VectorXd u(6); u.setZero() ; // the constrained velocity in generalized coordinate
    VectorXd u_tilde(6); // the unconstrained velocity
    u_tilde.block<3,1>(0,0) = m_objects[0]->xdot_tilde;
    u_tilde.block<3,1>(3,0) = m_objects[0]->omega_tilde;
    Vector3d &x = m_objects[0]->x;
    Quaterniond &q = m_objects[0]->q;

    // for plotting
    // Define matrices to store the values for each iteration
    MatrixXd lambda_n_values(num_contacts, m_settings.newtonIterations);
    MatrixXd lambda_f_values (2*num_contacts,m_settings.newtonIterations);
    MatrixXd delta_lambda_n_values(num_contacts, m_settings.newtonIterations);
    MatrixXd delta_lambda_f_values(2*num_contacts,m_settings.newtonIterations);
    MatrixXd u_values(u.size(), m_settings.newtonIterations);
    MatrixXd delta_u_values(u.size(), m_settings.newtonIterations);
    MatrixXd C_n_values(num_contacts,m_settings.newtonIterations);
    MatrixXd W_f_values(num_contacts,m_settings.newtonIterations);
    MatrixXd phi_n_values(num_contacts,m_settings.newtonIterations);



    ////////////////////////////////////////////////////////////////////////
    ////////////////// START THE NEWTON ITERATIONS ////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int nidx = 0; nidx < m_settings.newtonIterations; nidx++)
    {
        cout << "Newton It = " << nidx << endl;
        MatrixXd M_tilde = m_objects[0]->generalized_mass();
        MatrixXd M_inv_tilde = M_tilde.inverse();

        for (int cidx = 0; cidx < num_contacts; cidx++) {

            // preconditioning values
            double r_n = 1; double r_f = 1;

            cout << "Contact " << cidx << endl;

            // the collision
            auto &contact = contacts[cidx];

            // the contact basis
            MatrixXd n_t_b = contact->n_t_b;

            // the normal
            Vector3d n = n_t_b.col(0);
            MatrixXd t_b = n_t_b.block<3,2>(0,1);

            // object indices
            int objectA_idx = contact->objectA->m_idx;
            int objectB_idx = contact->objectB->m_idx;

            //
            double lambda_normal = lambda(cidx);
            Vector2d lambda_friction = lambda.block<2,1>(num_contacts + 2*cidx,0);

            cout << "lambda_normal = " << lambda_normal << endl;
            cout << "lambda_friction = " << lambda_friction << endl;

            // fill the normal constraint vector and normal jacobian
            if (objectA_idx == 0 )
            {

                // normal constraints
                double C_n = n.dot( ( x + q.toRotationMatrix()*contact->contact_wrt_objectA ) -
                        ( m_objects[1]->x + m_objects[1]->q.toRotationMatrix()*contact->contact_wrt_objectB ) ) - 1e-3;
                cout << "C_n (0) = " << C_n << " , n = " << n.transpose() << " , t = " << t_b.col(0).transpose() << " , b = " << t_b.col(1).transpose() << endl;

                if (C_n > 0)
                {
                    contact->m_is_active = false;
                }
                h(cidx) = fischer(C_n,lambda_normal) / dt;
                phi_n_values(cidx,nidx) = fischer(C_n,lambda_normal) / dt;
                double pfpa = fischer_partial_a(C_n, r_n * lambda_normal);
                double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal);
                cout << "pfpa = " << pfpa << endl;
                cout << "pfpb = " << pfpb << endl;
                C(cidx, cidx) = pfpb / (dt*dt) ;
                J.block<1, 3>(cidx, 0) = pfpa*n;
                J.block<1, 3>(cidx, 3) = -pfpa*n.cross(contact->contact_wrt_objectA);
                cout << "ContactwrtA = " << contact->contact_wrt_objectA.transpose() << endl;


                // friction constraints
                Matrix<double,6,2> D; D.block<3,2>(0,0) = t_b;
                D.block<3,1>(3,0) = -t_b.block<3,1>(0,0).transpose()*vec2skew(contact->contact_wrt_objectA);
                D.block<3,1>(3,1) = -t_b.block<3,1>(0,1).transpose()*vec2skew(contact->contact_wrt_objectA);
                J.block<2, 6>(num_contacts + 2 * cidx, 0) = D.transpose();
                double W_f = W_fischer(D,u, lambda_friction, r_f, 0.5, lambda_normal);
                h.block<2,1>(num_contacts + 2 * cidx,0) = D.transpose() * u + W_f * lambda_friction;
                C.block<2,2>(num_contacts + 2 * cidx,num_contacts + 2 * cidx) = W_f*Matrix2d::Identity() / dt ;

                C_n_values(cidx,nidx) = C_n;
                W_f_values(cidx,nidx) = W_f;
            }
            else if (objectA_idx == 1)
            {
                // normal constraints
                double C_n = n.dot(( m_objects[1]->x + m_objects[1]->q.toRotationMatrix()*contact->contact_wrt_objectA ) -
                        ( x + q.toRotationMatrix()*contact->contact_wrt_objectB ) ) - 1e-3;
                cout << "C_n (1) = " << C_n << " , n = " << n.transpose() << " , t = " << t_b.col(0).transpose() << " , b = " << t_b.col(1).transpose() << endl;

                if (C_n > 0)
                {
                    contact->m_is_active = false;
                }

                h(cidx) = fischer(C_n,lambda_normal) / dt;
                phi_n_values(cidx,nidx) = fischer(C_n,lambda_normal) / dt;
                double pfpa = fischer_partial_a(C_n, r_n * lambda_normal);
                double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal);
                cout << "pfpa = " << pfpa << endl;
                cout << "pfpb = " << pfpb << endl;
                C(cidx, cidx) = pfpb / (dt*dt) ;
                J.block<1, 3>(cidx, 0) = -pfpa*n;
                J.block<1, 3>(cidx, 3) = pfpa*n.cross(contact->contact_wrt_objectB);
                cout << "ContactwrtB = " << contact->contact_wrt_objectB.transpose() << endl;

                // friction constraints
                Matrix<double,6,2> D; D.block<3,2>(0,0) = -t_b;
                D.block<3,1>(3,0) = t_b.block<3,1>(0,0).transpose()*vec2skew(contact->contact_wrt_objectB);
                D.block<3,1>(3,1) = t_b.block<3,1>(0,1).transpose()*vec2skew(contact->contact_wrt_objectB);
                J.block<2, 6>( num_contacts + 2 * cidx, 0) = D.transpose();
                double W_f = W_fischer(D,u, lambda_friction, r_f, 0.5, lambda_normal);
                h.block<2,1>(num_contacts + 2*cidx,0) = D.transpose() * u + W_f * lambda_friction;
                C.block<2, 2>(2 * cidx, 2 * cidx) = W_f*Matrix2d::Identity() / dt ;

                C_n_values(cidx,nidx) = C_n;
                W_f_values(cidx,nidx) = W_f;
            }
        }

        MatrixXd H = M_tilde;
        VectorXd g = M_tilde * (u - u_tilde) - dt * J.transpose() * lambda;

        auto H_inv = H.inverse();
        MatrixXd A = J * H_inv * J.transpose() + C;
        A += 1e-6*MatrixXd::Identity(A.rows(),A.cols());

        VectorXd b = (J * H_inv * g - h) / dt;

        cout << "J = " << J << endl;

        // zero out any inactive contacts
        for (int cidx = 0; cidx < num_contacts; cidx++)
        {
            if (m_collisionDetector->m_contacts[cidx]->is_active() == 0 )
            {
                b(num_contacts + 2*cidx) = 0;
                b(num_contacts + 2*cidx + 1) = 0;
            }
        }

        // solve this system with Gauss-Seidel
        VectorXd delta_lambda = A.inverse()*b;

        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
        // perform a line search
        double t = 0.5;

        // this function computes the system residual
        auto r = computeResidual(x,q,u,u_tilde,lambda,M_tilde,contacts,dt);
        cout << "r = " << r.transpose() << endl;

        // update solution
        lambda += t*delta_lambda;
        u += t*delta_u;
        x += delta_u.block<3,1>(0,0)*dt;
        q = sumQuaternions(q,
                       multiplyQuaternionScalar(dt * 0.5 , angularVelocityToQuaternion(delta_u.block<3,1>(3,0)) * q));
        q.normalize();
        m_objects[0]->update_inertia_matrix();

        // Store the values for this iteration
        lambda_n_values.col(nidx) = lambda.head(num_contacts);
        lambda_f_values.col(nidx) = lambda.tail(2*num_contacts);
        delta_lambda_n_values.col(nidx) = delta_lambda.head(num_contacts);
        delta_lambda_f_values.col(nidx) = delta_lambda.tail(2*num_contacts);
        u_values.col(nidx) = u;
        delta_u_values.col(nidx) = delta_u;

        cout << "Delta Lambda : " << endl << delta_lambda.transpose() << endl;
        cout << "Delta U : " << endl << delta_u.transpose() << endl;
//        cout << "lambda : " << endl << lambda.transpose() << endl;
//        cout << "u : " << u.transpose() << endl;

    }
    m_objects[0]->x = x;
    m_objects[0]->q = q;
    m_objects[0]->xdot = u.block<3,1>(0,0);
    m_objects[0]->omega = u.block<3,1>(3,0);

    // Create a 2x2 grid of subplots
    plt::subplot(3, 3, 1);
    plt::title("Lambda N Values");
    for (int i = 0; i < lambda_n_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = lambda_n_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::named_plot(legend.str(), row_i_std);
    }
    plt::legend();

    plt::subplot(3, 3, 2);
    plt::title("Lambda F Values");
    for (int i = 0; i < lambda_f_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = lambda_f_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::named_plot(legend.str(), row_i_std);
    }
    plt::legend();

    plt::subplot(3, 3, 3);
    plt::title("Delta Lambda N Values");
    for (int i = 0; i < delta_lambda_n_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = delta_lambda_n_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 4);
    plt::title("Delta Lambda F Values");
    for (int i = 0; i < delta_lambda_f_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = delta_lambda_f_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 5);
    plt::title("U Values");
    for (int i = 0; i < u_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = u_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 6);
    plt::title("Delta U Values");
    for (int i = 0; i < delta_u_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = delta_u_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 7);
    plt::title("C_n Values");
    for (int i = 0; i < C_n_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = C_n_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 8);
    plt::title("W_f Values");
    for (int i = 0; i < W_f_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = W_f_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::subplot(3, 3, 9);
    plt::title("Phi_n Values");
    for (int i = 0; i < phi_n_values.rows(); ++i) {
        std::stringstream legend;
        legend << "Row " << i;

        VectorXd row_i = phi_n_values.row(i);
        std::vector<double> row_i_std(row_i.data(), row_i.data() + row_i.size());

        plt::plot(row_i_std, {{"label", legend.str()}});
    }
    plt::legend();

    plt::show();


    cout << "End Newton It " << endl;
}