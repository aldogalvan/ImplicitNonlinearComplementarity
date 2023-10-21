
#include "PenaltyDemo.hpp"


inline double W_fischer(const MatrixXd& D, const VectorXd& qdot, const VectorXd& lambda_f,
                        const double& r, const double& mu, const double& lambda_n)
{
    double ret = 1;
    if (lambda_n * mu > 0)
    {
        ret *= r*(sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)) - r*(mu*lambda_n - lambda_f.norm()))/
               ((D.transpose()*qdot).squaredNorm()+r*mu*lambda_n-sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)));
    }
    if (isnan(ret))
        ret = 0;

    return ret;
}

inline double fischer_partial_a(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    double ret;
    if (abs(a) < 1e-12 && abs(b) < 1e-12)
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
    if (abs(a) < 1e-12 && abs(b) < 1e-12)
    {
        ret = 1;
    }
    else
    {
        ret = 1 - (b/sqrt(pow(a,2) + pow(b,2)));
    }
    return ret;
}

inline double fischer(double C_n, double lambda_n) {
    return  C_n + lambda_n - sqrt(C_n * C_n + lambda_n * lambda_n) ;
}

VectorXd PenaltyDemo::computeResidual(const VectorXd &x, const VectorXd &u,
                                    const VectorXd &delta_u, const VectorXd &u_tilde,
                                    const VectorXd& lambda, const VectorXd& delta_lambda, const MatrixXd &M_tilde,
                                    const vector<Contact*>& contacts, double dt) {
    int num_contacts = contacts.size();
    VectorXd r(3 + lambda.rows());
    MatrixXd J_n(num_contacts,3);
    MatrixXd J_f(2*num_contacts,3);
    VectorXd lambda_n = lambda.head(num_contacts);
//    VectorXd lambda_f = lambda.tail(2*num_contacts);

    for (int cidx = 0 ; cidx < num_contacts; cidx ++)
    {
        Contact* contact = contacts[cidx];

        // normal contact constraints
        double lambda_normal = lambda_n(cidx);
        Vector3d n = contact->n_t_b.col(0);
        double C_n;
        if (contact->objectIdxB == -1)
        {
            C_n = -n.dot((contact->contact_wrt_objectB + block->getLocalPos().eigen()) - (contact->contact_wrt_objectA + x));
        }
        else if (contact->objectIdxB == 0)
        {
            C_n = -n.dot((contact->contact_wrt_objectB + x) - (contact->contact_wrt_objectA + block->getLocalPos().eigen()));
        }

        VectorXd pphin_pq = fischer_partial_a(C_n,lambda_normal)*n;
        J_n.row(cidx) = pphin_pq.transpose();
        r(3 + cidx) = fischer(C_n,lambda_normal);

        // friction contact constraints
//                double r = 1; double mu = 0.5;
//                Vector2d lambda_friction = lambda.block<2,1>(3*cidx+1,0);
//                Matrix<double,3,2> D = contact->n_t_b.block<3,2>(0,1);
//                MatrixXd pphif_pqdot = D;
//                double W_f = W_fischer(D,u,lambda_friction,r,mu,lambda_normal);
//                MatrixXd pphif_plambda = W_f*Matrix2d::Identity();
//                MatrixXd phi_friction = D.transpose()*u + W_f*lambda_friction;
//                J.block<2,3>(3*cidx+1,0) = pphif_pqdot.transpose();
//                C.block<2,2>(3*cidx + 1, 3*cidx + 1) = pphif_plambda / dt;
    }
    r.block<3,1>(0,0) = M_tilde*(u - u_tilde)/dt + J_n.transpose()*lambda_n;

    return r;
}

void PenaltyDemo::backtrackingLineSearch(double &t, double &err,
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

void PenaltyDemo::initialize()
{
    // create the meshes
//    peg = new GodObject(m_devicePtr,0,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/RigidBodyDemo/cube.obj");
//    peg->scaleMesh(2);
    peg = new GodObject(m_devicePtr,0,"/home/aldo/ImplicitNonlinearComplementarity/resources/PegInHole/peg.obj");
    peg->scaleMesh(0.025);
    peg->scaleMesh(0.99);
    peg->getMesh(0)->m_material->setBlue(); peg->getMesh(0)->m_material->setShininess(0.5);
    m_world->addChild(peg);
    peg->compute_inertia_matrix();
    peg->setWireMode(true);

    cout << "Peg Vertices = " << peg->getNumVertices() << endl;

    // create the block
//    block = new RigidObject(0,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/RigidBodyDemo/cube.obj");
//    block->scaleMesh(2);
    block = new RigidObject(0,"/home/aldo/ImplicitNonlinearComplementarity/resources/PegInHole/block.obj");
    block->scaleMesh(.025);
    block->getMesh(0)->m_material->setRed(); block->getMesh(0)->m_material->setShininess(0.5);
    m_world->addChild(block);
    block->compute_inertia_matrix();
    block->set_local_pos(Vector3d(0,0,-0.45));
    block->computeAllNormals();
    block->setWireMode(true);

    cout << "Block Vertices = " << block->getNumVertices() << endl;
}

void PenaltyDemo::updateHaptics(Vector3d& f)
{
    peg->updateFromDevice();
}

void PenaltyDemo::step(double dt)
{

    // update the variables of the complementarity problem
    peg->updateComplementarityProblem(dt);

    // update the block
    block->x_tilde = block->x;

    VectorXd peg_start = peg->get_configuration();
    VectorXd peg_end = peg->get_configuration_unconstrained();
    const MatrixXd& peg_vertices = peg->vertices();
    const MatrixXi& peg_triangles = peg->triangles();
    VectorXd block_start = block->get_configuration();
    VectorXd block_end = block->get_configuration_unconstrained();
    const MatrixXd& block_vertices = block->vertices();
    const MatrixXi& block_triangles = block->triangles();

    vector<Contact*> contacts;
    vector<cShapeLine*> colvis;

    if (CollisionDetector::findCollisionsRigidRigid(peg_start,peg_end,peg_vertices, peg_triangles,
                                                    block_start, block_end, block_vertices, block_triangles, contacts,colvis))
    {

        for (int cidx = 0 ; cidx < colvis.size() ; cidx ++)
        {    m_world->addChild(colvis[cidx]);}

        int ndof = 3; // the degrees of freedom of the god object
        int num_contacts = contacts.size();

        cout << "Num Contacts = " << num_contacts << endl;

        // God Object
        Vector3d x = peg->x_tilde;
        Vector3d u = peg->xdot_tilde;
        MatrixXd M = peg->M();
        const Vector3d& u_tilde = peg->xdot_tilde;

        // Block
        Vector3d& x_block = block->x;

        // lagrange multipliers
        VectorXd lambda =  VectorXd::Zero(num_contacts + 2*num_contacts);
        VectorXd delta_lambda = VectorXd::Zero(num_contacts + 2*num_contacts);
        VectorXd h = VectorXd::Zero(num_contacts + 2*num_contacts);
        MatrixXd C = MatrixXd::Zero(num_contacts + 2*num_contacts,num_contacts + 2*num_contacts);
        VectorXi is_active = VectorXi::Zero(num_contacts + 2*num_contacts);

        // iterates over contacts to create the contact jacobian
        J.resize(num_contacts + 2*num_contacts,ndof);

        // compute the current residual

        for (int newtonit = 0; newtonit < numNewtonIt; newtonit++)
        {
            cout << "Newton It : " << newtonit << endl;
            cout << "x = " << x.transpose() << endl;
            cout << "u = " << u.transpose() << endl;
            cout << "u_tilde = " << u_tilde.transpose() << endl;
            for (int cidx = 0; cidx < num_contacts; cidx++)
            {

                // the contact of interest
                Contact* contact = contacts[cidx];

                // normal contact constraints
                double lambda_normal = lambda(cidx);

                // friction
                Vector2d lambda_friction = lambda.block<2,1>(num_contacts + 2 * cidx, 0);

                // the contact basis
                MatrixXd n_t_b = contact->n_t_b;

                // the normal
                Vector3d n = n_t_b.col(0);
                MatrixXd D = n_t_b.block<3,2>(0,1);

                // object indices
                int objectA_idx = contact->objectIdxA;
                int objectB_idx = contact->objectIdxB;

                // preconditioning parameter
                double r_n = 1; double r_f = 1;

                // fill the normal constraint vector and normal jacobian
                if (objectA_idx == 0 )
                {
                    // normal constraints
                    double C_n = n.dot(( x + contact->contact_wrt_objectA ) - (block->x + contact->contact_wrt_objectB) ) - 1e-6;
                    cout << "C_n = " << C_n << endl;
                    h(cidx) = fischer(C_n,lambda_normal) / dt;
                    double pfpa = fischer_partial_a(C_n, r_n * lambda_normal);
                    double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal);
                    C(cidx, cidx) = pfpb / (dt*dt) ;
                    J.block<1, 3>(cidx, 0) = pfpa*n;

                    // friction constraints
                    J.block<2, 3>(num_contacts + 2 * cidx, 0) = D.transpose();
                    double W_f = W_fischer(D,u, lambda_friction, r_f, 0.5, lambda_normal);
                    h.block<2,1>(num_contacts + 2*cidx,0) = D.transpose() * u + W_f * lambda_friction;
                    C.block<2, 2>(num_contacts + 2 * cidx, num_contacts + 2 * cidx) = W_f*Matrix2d::Identity() / dt ;
                }
                else if (objectA_idx == 1)
                {
                    // normal constraints
                    double C_n = n.dot(( block->x + contact->contact_wrt_objectA ) - ( x + contact->contact_wrt_objectB ) );
                    cout << "C_n = " << C_n << endl;
                    h(cidx) = fischer(C_n,lambda_normal) / dt;
                    double pfpa = fischer_partial_a(C_n, r_n * lambda_normal);
                    double pfpb =  fischer_partial_b(C_n, r_n * lambda_normal);
                    C(cidx, cidx) = pfpb / (dt*dt) ;
                    J.block<1, 3>(cidx, 0) = -pfpa*n;

                    // friction constraints
                    J.block<2, 3>( num_contacts + 2 * cidx, 0) = -D.transpose();
                    double W_f = W_fischer(D,u, lambda_friction, r_f, 0.5, lambda_normal);
                    h.block<2 , 1>(num_contacts + 2*cidx,0) = D.transpose() * u + W_f * lambda_friction;
                    C.block<2 , 2>(num_contacts + 2 * cidx, num_contacts + 2 * cidx) = W_f*Matrix2d::Identity() / dt ;
                }
            }

//            cout << "J = " << J << endl;

            // Assemble system matrices
            MatrixXd H = M; MatrixXd H_inv = H.inverse();
            VectorXd g = M * (u - u_tilde) - dt * J.transpose() * lambda;
            MatrixXd A = J * H_inv * J.transpose() + C; A += 1e-6*MatrixXd::Identity(A.rows(),A.cols());
            VectorXd b = (J * H_inv * g - h) / dt;

            for (int cidx = 0 ; cidx < num_contacts; cidx++)
                1;

            // solve the system
            Eigen::FullPivLU<Eigen::MatrixXd> lu_decomp(A); // Perform LU decomposition
            delta_lambda = lu_decomp.solve(b);
            VectorXd delta_u = (H_inv * (J.transpose() * (delta_lambda * dt) - g));
            VectorXd impulse = H_inv * (J.transpose() * (delta_lambda * dt));
            cout << "impulse = " << impulse.transpose() << endl;
            cout << "delta_u = " << delta_u.transpose() << endl;

            // perform a backtracking line search
            double t = 0.75;

            for (int i = 0 ; i < delta_lambda.size(); i++)
            {
                if (delta_lambda[i] < 0 )
                    delta_lambda[i] = 0;
            }

            // the residual of this function
            // VectorXd r = computeResidual(x, u, u_tilde, lambda_normal, M, colInfo, dt);

            // the backtracking line search
//            backtrackingLineSearch(t, err, x , u , u_tilde, delta_u, lambda_normal , delta_lambda,
//                                   M, colInfo, dt, alpha, beta, backtrackingIt );

            cout << "delta_lambda = " << delta_lambda.transpose() << endl;
            lambda += delta_lambda*t;
            u += delta_u*t;
            cout << "delta_x = " << delta_u.transpose()*dt << endl;
            x += delta_u*dt;

            auto res = computeResidual(x,u,delta_u,u_tilde,lambda,delta_lambda,M,contacts,dt);
            cout << "Residual = " << res.transpose() << endl;

        }
        peg->x = x;
        peg->xdot = u;

        for (int cidx = 0 ; cidx < contacts.size(); cidx++)
        {
            delete contacts[cidx];
        }

        cout << "END" << endl;
    }
    else
    {
        peg->x = peg->x_tilde;
        C.resize(0,0); J.resize(0,0);
    }

}

void PenaltyDemo::updateGraphics()
{
    peg->update_mesh_position();
    block->update_mesh_position();
    vis->setLocalPos(peg->x_d);
}
