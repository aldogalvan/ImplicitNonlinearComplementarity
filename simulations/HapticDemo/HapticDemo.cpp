#include "HapticDemo.h"

inline MatrixXd tangent_basis(const Vector3d& normal)
{
    MatrixXd D(3, 2);

    // Normalize the normal vector
    Vector3d normalizedN = normal.normalized();

    // random
    Vector3d random_vector = Vector3d::Random();
    random_vector.normalize();

    // first random vector
    Vector3d V1 = normalizedN.cross(random_vector).normalized();
    assert(abs(V1.dot(normalizedN)) < 1e-12);

    // Compute the second direction vector using cross product
    Vector3d V2 = normalizedN.cross(V1).normalized();

    assert(abs(V2.dot(normalizedN)) < 1e-12);

    // Store the direction vectors in the matrix D
    D.col(0) = V1;
    D.col(1) = V2;

    return D;
}

inline double W_fischer(const MatrixXd& D, const VectorXd& qdot, const VectorXd& lambda_f,
                        const double& r, const double& mu, const double& lambda_n)
{
    double ret = 1;
    if (lambda_n * mu > 0)
    {
        ret *= r*(sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)) - r*(mu*lambda_n - lambda_f.norm()))/
               ((D.transpose()*qdot).squaredNorm()+r*mu*lambda_n-sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)));
        if (isnan(ret))
            ret = 1;
    }
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


inline VectorXd computeResidual(const Vector3d &x, const Vector3d &u, const Vector3d &u_tilde,
                                VectorXd lambda_n, const Matrix3d &M_tilde, const ColInfo* colInfo , double dt) {
    VectorXd r(M_tilde.rows() + 1, 1);
    int start_idx = M_tilde.rows();

    r.head(M_tilde.rows()) = M_tilde * (u - u_tilde) / dt - colInfo->normal * lambda_n;
    double C_n = colInfo->normal.dot(x - colInfo->contactPt);
    r(3) = fischer(C_n, lambda_n(0));

    return r;
}

static inline void backtrackingLineSearch(double &t, double &err,
                                          const VectorXd &x, const VectorXd &u,
                                          const VectorXd &u_tilde, const VectorXd &delta_u,
                                          VectorXd lambda_n, VectorXd delta_lambda_n,
                                          const MatrixXd &M_tilde, const ColInfo* colInfo, double dt, double alpha, double beta, int max_iter) {
    VectorXd grad_err = delta_lambda_n;
    double err_k;

    for (int k = 0; k < max_iter; ++k) {
        VectorXd lambda_n_k = lambda_n + t * delta_lambda_n;
        VectorXd u_k = u + t * delta_u;
        VectorXd x_k = x + u * dt;
        VectorXd r = computeResidual(x_k, u_k, u_tilde, lambda_n_k, M_tilde, colInfo, dt);
        err_k = 0.5 * r.dot(r);

        // Perform Armijo condition to check if we have sufficient decrease
        if (err_k <= err + beta * t * grad_err.dot(delta_lambda_n)) {
            break;
        }

        t = alpha * t;
    }

    err = err_k;
}

inline bool computeCollision(Plane* plane, GodObject* godObj, ColInfo* colInfo)
{
    // Compute the distance between the sphere center and the plane
    Eigen::Vector3d sphereCenter = godObj->x_tilde;
    double distance = (sphereCenter - plane->point).dot(plane->normal);

    // Check for collision
    if (distance <= godObj->r)
    {
        colInfo->normal = plane->normal;
        colInfo->contactPt = Vector3d(sphereCenter(0),sphereCenter(1),0) + plane->point;
        colInfo->depth = godObj->r - distance;
        return true;
    }
    else
    {
        return false;
    }
}

void HapticDemo::updateHaptics(Vector3d& f, double dt) {

    // update from the device
    godObject->updateFromDevice();

    int ndof = 3;
    int numContacts = 1;

    // legibility
    Vector3d x = godObject->x_d;
    Vector3d u = Vector3d(0,0,0);
    Matrix3d &M = godObject->M;
    Matrix3d B = Matrix3d::Zero();
    Matrix3d K = Matrix3d::Zero();
    const Vector3d &u_tilde = Vector3d(0,0,0);

    // lagrange multipliers
    VectorXd lambda = VectorXd::Zero(numContacts + 2 * numContacts);
    VectorXd delta_lambda = VectorXd::Zero(numContacts + 2 * numContacts);
    MatrixXd J_normal;
    MatrixXd J_friction;
    MatrixXd D;
    Vector2d phi_friction = Vector2d::Zero();
    double phi_normal = 0;
    bool is_active = 0;
    double S = 0; MatrixXd W = Matrix2d::Zero();
    if (colInfo != NULL)
        D = tangent_basis(colInfo->normal);

    for (int newtonit = 0; newtonit < numNewtonIt; newtonit++) {
        cout << "Newton It = " << newtonit << endl;
        cout << "x = " << x.transpose() << endl;
//            cout << "u = " << u.transpose() << endl;
//            cout << "u_tilde = " << u_tilde.transpose() << endl;

        if (colInfo != NULL)
        {
            // normal
            double lambda_normal = lambda(0);
            double C_n = colInfo->normal.dot((x - Vector3d(0, 0, 1) * godObject->r) - colInfo->contactPt) - contactTol;
//            cout << "C_n = " << C_n << endl;
            if (C_n < -contactTol && delta_lambda(0) > 0) is_active = 1;
            VectorXd pphin_pq = fischer_partial_a(C_n, lambda_normal) * colInfo->normal;
            double pphin_plambda = fischer_partial_b(C_n, lambda_normal);
            phi_normal = fischer(C_n, lambda_normal);
            J_normal = pphin_pq.transpose();
            S = pphin_plambda;
            phi_normal /= dt;

            // friction constraints
            double r = 1;
            double mu = 0.0;
            VectorXd lambda_friction = lambda.tail(2);
            MatrixXd pphif_pqdot = D;
            double W_f = W_fischer(D, u, lambda_friction, r, mu, lambda_normal);
            MatrixXd pphif_plambda = W_f * Matrix2d::Identity();
            phi_friction = D.transpose() * u + W_f * lambda_friction;
            J_friction = pphif_pqdot.transpose();
            W = pphif_plambda;
        }

        // Assemble system matrices
        MatrixXd H = M;
        MatrixXd H_inv = H.inverse();
        MatrixXd J = Matrix3d::Zero();
        J.row(0) = J_normal;
        J.block<2, 3>(1, 0) = J_friction;
        MatrixXd C = MatrixXd::Zero(3, 3);
        C(0, 0) = S / (dt * dt);
        C.block<2, 2>(1, 1) = W / dt;
        VectorXd g = M * (u - u_tilde) - dt * J.transpose() * lambda;
        VectorXd h(3);
        h(0) = phi_normal;
        h.tail(2) = phi_friction;
        MatrixXd A = J * H_inv * J.transpose() + C;
        A += 1e-6 * MatrixXd::Identity(A.rows(), A.cols());
        VectorXd b = (J * H_inv * g - h) / dt;

        if (!is_active) {
//                A.block<2,3>(1,0).setZero();
            b.tail(2).setZero();
        }

        // solve the system
        FullPivLU<MatrixXd> lu_decomp(A); // Perform LU decomposition
        delta_lambda = lu_decomp.solve(b);
//            cout << "delta_lambda = " << delta_lambda.transpose() << endl;
        VectorXd delta_u = (H_inv * (J.transpose() * (delta_lambda * dt) - g));

        // perform a backtracking line search
        double t = 1;

        // the residual of this function
        // VectorXd r = computeResidual(x, u, u_tilde, lambda_normal, M, colInfo, dt);

        // the backtracking line search
//            backtrackingLineSearch(t, err, x , u , u_tilde, delta_u, lambda_normal , delta_lambda,
//                                   M, colInfo, dt, alpha, beta, backtrackingIt );

//            cout << "delta_u" << delta_u << endl;
        lambda += delta_lambda * t;
        u += delta_u * t;
        x += delta_u * dt;

    }

    // PID haptic control
    f = godObject->computeForce();
}

void HapticDemo::step(double dt)
{
    colInfo = new ColInfo();

    // update the variables of the complementarity problem
    godObject->updateComplementarityProblem(dt);

    if (computeCollision(plane,godObject, colInfo))
    {
        int ndof = 3;
        int numContacts = 1;

        // legibility
        Vector3d x = godObject->x;
        Vector3d u = godObject->u;
        Matrix3d& M = godObject->M;
        const Vector3d& u_tilde = godObject->u_tilde;

        // lagrange multipliers
        VectorXd lambda =  VectorXd::Zero(numContacts + 2*numContacts);
        VectorXd delta_lambda = VectorXd::Zero(numContacts + 2*numContacts);
        MatrixXd D = tangent_basis(colInfo->normal);

        for (int newtonit = 0; newtonit < numNewtonIt; newtonit++)
        {
            cout << "Newton It = " << newtonit << endl;
            cout << "x = " << x.transpose() << endl;
//            cout << "u = " << u.transpose() << endl;
//            cout << "u_tilde = " << u_tilde.transpose() << endl;

            // normal
            double lambda_normal = lambda(0);
            double C_n = colInfo->normal.dot((x - Vector3d(0,0,1)*godObject->r) - colInfo->contactPt) - contactTol;
//            cout << "C_n = " << C_n << endl;
            bool is_active = 0; if (C_n < -contactTol && delta_lambda(0) > 0) is_active = 1;
            VectorXd pphin_pq = fischer_partial_a(C_n,lambda_normal)*colInfo->normal;
            double pphin_plambda = fischer_partial_b(C_n,lambda_normal);
            double phi_normal = fischer(C_n,lambda_normal);
            MatrixXd J_normal = pphin_pq.transpose();
            double S = pphin_plambda;
            phi_normal /= dt;

            // friction constraints
            double r = 1; double mu = 0.0;
            VectorXd lambda_friction = lambda.tail(2);
            MatrixXd pphif_pqdot = D;
            double W_f = W_fischer(D,u,lambda_friction,r,mu,lambda_normal);
            MatrixXd pphif_plambda = W_f*Matrix2d::Identity();
            MatrixXd phi_friction = D.transpose()*u + W_f*lambda_friction;
            MatrixXd J_friction = pphif_pqdot.transpose();
            MatrixXd W = pphif_plambda;


            // Assemble system matrices
            MatrixXd H = M; MatrixXd H_inv = H.inverse();
            if (J.rows() != J_normal.rows() + J_friction.rows())
                J.resize(J_normal.rows() + J_friction.rows(),ndof);
            J.row(0) = J_normal; J.block<2,3>(1,0) = J_friction;

            MatrixXd C = MatrixXd::Zero(3,3); C(0,0) = S / (dt * dt); C.block<2,2>(1,1) = W / dt;
            VectorXd g = M * (u - u_tilde) - dt * J.transpose() * lambda;
            VectorXd h(3); h(0) = phi_normal; h.tail(2) = phi_friction;
            C_delassus = J * H_inv * J.transpose();
            MatrixXd A = C_delassus + C; A += 1e-6*MatrixXd::Identity(A.rows(),A.cols());
            VectorXd b = (J * H_inv * g - h) / dt;

            if (!is_active)
            {
//                A.block<2,3>(1,0).setZero();
                b.tail(2).setZero();
            }

            // solve the system
            FullPivLU<MatrixXd> lu_decomp(A); // Perform LU decomposition
            delta_lambda = lu_decomp.solve(b);
//            cout << "delta_lambda = " << delta_lambda.transpose() << endl;
            VectorXd delta_u = (H_inv * (J.transpose() * (delta_lambda * dt) - g));

            // perform a backtracking line search
            double t = 1;

            // the residual of this function
            // VectorXd r = computeResidual(x, u, u_tilde, lambda_normal, M, colInfo, dt);

            // the backtracking line search
//            backtrackingLineSearch(t, err, x , u , u_tilde, delta_u, lambda_normal , delta_lambda,
//                                   M, colInfo, dt, alpha, beta, backtrackingIt );

//            cout << "delta_u" << delta_u << endl;
            lambda += delta_lambda*t;
            u += delta_u*t;
            x += delta_u*dt;

        }
        godObject->x = x;
        godObject->u = u;
        cout << "END" << endl;
    }
    else
    {
        godObject->x = godObject->x_d;
        J.resize(0,0);
    }

    delete colInfo;
}

void HapticDemo::updateGraphics()
{
    godObject->updateVisualizer();
}
