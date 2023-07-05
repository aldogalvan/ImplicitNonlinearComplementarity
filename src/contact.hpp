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

    // time of collision
    double t;
    // collision normal
    // defined with respect to bodyB
    // the 'triangle owner'
    Vector3d n;
    // collision depth
    double depth;
    // collision body indices
    int bodyIdxA, bodyIdxB;
    // the global position of the contact
    Vector3d contact_pt;
    // the contact point for object A
    // defined in rigid body frame
    Vector3d contact_wrt_bodyA;
    // the contact pont for object B
    // defined in rigid body frame
    Vector3d contact_wrt_bodyB;
    // the pointer to object A
    RigidObject* objectA;
    RigidObject* objectB;
};


#endif //IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP
