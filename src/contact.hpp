#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP

struct Contact
{
    // time of collision
    double t;
    // collision normal
    // defined with respect to bodyB
    // the 'triangle owner'
    Vector3d normal;
    // collision depth
    double depth;
    // collision body indices
    int bodyIdxA, bodyIdxB;
    // the position of the ocntact
    Vector3d contact_pt;
    // the contact point for object A
    // defined in rigid body frame
    Vector3d contact_wrt_bodyA;
    // the contact pont for object B
    // defined in rigid body frame
    Vector3d contact_wrt_bodyB;
    // the pointer to object A
    Object* objectA;
    Object* objectB;
};

struct ContactDeformable
{

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_CONTACT_HPP
