#include "contact.hpp"

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
    }
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

void Contact::compute_contact_basis()
{
    n_t_b.resize(3,3);
    n_t_b.col(0) = n;
    auto t_b = tangent_basis(n);
    n_t_b.col(1) = t_b.col(0);
    n_t_b.col(2) = t_b.col(1);
}

void Contact::compute_constraints(double r_n, double r_f, double lambda_n, VectorXd lambda_f)
{

    // collision values
    Vector3d normal = this->n;
    Vector3d pt = this->contact_pt;

    // body A configuration
    Vector3d pos_A = objectA->x;
    Quaterniond rot_A = objectA->q;
    Vector3d xdot_A = objectA->xdot;
    Vector3d omega_A = objectA->omega;
    rot_A.normalize();

    // body B configuration
    Vector3d pos_B = objectB->x;
    Quaterniond rot_B = objectB->q;
    rot_B.normalize();
    Vector3d xdot_B = objectB->xdot;
    Vector3d omega_B = objectB->omega;

    // the contact pt wrt A
    VectorXd contact_ptA = rot_A*contact_wrt_objectA + pos_A;
    // the contact pt wrt A
    VectorXd contact_ptB = rot_B*contact_wrt_objectB + pos_B;

    // compute penetration depth (should already be done but w/e)
    double penetration_depth = (contact_ptA - contact_ptB).dot(normal);

    // if the penetration depth is > 0 constraint is inactive
    if (penetration_depth > 0) {
        m_is_active = 0;
    } else {
        m_is_active = 1;
    }

    // the kinematic map matrices
    MatrixXd G_A = objectA->kinematic_map_G();
    MatrixXd G_B = objectB->kinematic_map_G();

    ///////////////////////////////////////////////////////////////////////
    ////// Assemble the Normal Constraint Matrices and Jacobians /////////
    //////////////////////////////////////////////////////////////////////

    // the constraint value
    double C_n = penetration_depth;

//    cout << "C_n = " << C_n << endl;
//    cout << "normal = " << normal.transpose() << endl;

    //
    // gradients of the constraint for normal
    VectorXd grad_C_n_A(6); grad_C_n_A.head(3) = normal; grad_C_n_A.tail(3) = vec2skew(-normal).transpose()*contact_wrt_objectA;
    VectorXd grad_C_n_B(6); grad_C_n_B.head(3) = -normal; grad_C_n_B.tail(3) = vec2skew(normal).transpose()*contact_wrt_objectB;

    // fill the normal constraint value
    phi_n = fischer(C_n, r_n * lambda_n);

    double d_fischer_d_a = fischer_partial_a(C_n, r_n * lambda_n);

    J_n_A = (d_fischer_d_a * grad_C_n_A).transpose();
//    cout << "J_n_A = " << J_n_A << endl;
    assert(!J_n_A.hasNaN());
    J_n_B = (d_fischer_d_a * grad_C_n_B).transpose();
//    cout << "J_n_B = " << J_n_B << endl;
    assert(!J_n_B.hasNaN());
    phi_n_wrt_lambda_n = fischer_partial_b(C_n, r_n * lambda_n) * r_n;
    assert(!isnan(phi_n_wrt_lambda_n));

    ///////////////////////////////////////////////////////////////////////
    ////// Assemble the Friction Constraint Matrices and Jacobians ///////
    //////////////////////////////////////////////////////////////////////

    // the relative velocity
    Vector3d qdot = (xdot_A + omega_A.cross(contact_wrt_objectA)) -
            (xdot_B + omega_B.cross(contact_wrt_objectA));

    MatrixXd D = n_t_b.block<3,2>(0,1);
    // fill the friction vector
    double W_f = W_fischer(D,qdot, lambda_f, r_f, mu, lambda_n);
    phi_f = D.transpose() * qdot+ W_f * lambda_f;
    auto phi_friction_wrt_qdot = D;
    phi_f_wrt_lambda_f = W_f * Matrix2d::Identity();
    Vector3d t = n_t_b.col(1); Vector3d b = n_t_b.col(2);
    MatrixXd D_n_A(6,2); MatrixXd D_n_B(6,2);
    D_n_A.block<3,1>(0,0) = t; D_n_A.block<3,1>(0,1) = b;
    D_n_A.block<3,1>(3,0)= vec2skew(-t).transpose()*(contact_wrt_objectA);
    D_n_A.block<3,1>(3,1)= vec2skew(-b).transpose()*(contact_wrt_objectA);
    D_n_B.block<3,1>(0,0) = -t; D_n_A.block<3,1>(0,1) = -b;
    D_n_B.block<3,1>(3,0)= vec2skew(t).transpose()*contact_wrt_objectB;
    D_n_B.block<3,1>(3,1)= vec2skew(b).transpose()*contact_wrt_objectB;
    J_f_A = D_n_A.transpose(); J_f_B = D_n_B.transpose();

}

MatrixXd Contact::Jacobian_normal_A()
{
    return J_n_A;
}

MatrixXd Contact::Jacobian_normal_B()
{
    return J_n_B;
}

MatrixXd Contact::Jacobian_friction_A()
{
    return J_f_A;
}

MatrixXd Contact::Jacobian_friction_B()
{
    return J_f_B;
}

double Contact::phi_normal()
{
    return phi_n;
}

VectorXd Contact::phi_friction()
{
    return phi_f;
}

double Contact::partial_phi_normal_wrt_lambda_normal()
{
    return phi_n_wrt_lambda_n;
}

MatrixXd Contact::partial_phi_friction_wrt_lambda_friction()
{
    return phi_f_wrt_lambda_f;
}

