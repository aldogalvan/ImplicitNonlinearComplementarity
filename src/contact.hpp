#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP

#include <Eigen/Dense>
#include "objects.hpp"

using namespace Eigen;
class Contact
{
public:

    Contact(RigidObject* _objectA, RigidObject* _objectB)
    {
        objectA = _objectA; objectB = _objectB;
        objectA->m_contacts.emplace_back(this); objectB->m_contacts.emplace_back(this);
    }

    Contact(){}
    ~Contact() = default;

    void compute_constraints(double r_n , double r_f,
                             double lambda_n, VectorXd lambda_f);
    void compute_contact_basis();

    MatrixXd Jacobian_normal_A();
    MatrixXd Jacobian_normal_B();
    MatrixXd Jacobian_friction_A();
    MatrixXd Jacobian_friction_B();
    double phi_normal();
    VectorXd phi_friction();
    double partial_phi_normal_wrt_lambda_normal();
    MatrixXd partial_phi_friction_wrt_lambda_friction();
    bool is_active(){return m_is_active;}

    double phi_n; VectorXd phi_f; // the constraint values for object A
    MatrixXd J_n_A; MatrixXd J_n_B; // the normal jacobians
    MatrixXd J_f_A; MatrixXd J_f_B; // the friction jacobian
    double phi_n_wrt_lambda_n;
    MatrixXd phi_f_wrt_lambda_f;
    bool m_is_active; // is contact active
    MatrixXd n_t_b; // the contact basis (per column)

    double t;     // time of collision
    Vector3d n;     // collision normal from B to A
    double depth;     // collision depth
    int objectIdxA, objectIdxB;     // collision body indices
    Vector3d contact_pt;     // the global position of the contact
    Vector3d contact_wrt_objectA;    // the contact point for object A
    Vector3d contact_wrt_objectB;     // the contact pont for object B
    RigidObject* objectA;     // the pointer to object A
    RigidObject* objectB;     // the pointer to object B
    double mu = 0.5;            // friciton constant
    double contactTol = 1e-6; // contact tolerance
};


#endif //IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP
