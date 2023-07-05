
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H

#include "chai3d.h"
#include <Eigen/Dense>
#include "../../src/objects.hpp"
#include "../../src/collision.hpp"
#include "../../src/helper.hpp"

struct Settings
{

};

class RigidBodyDemo
{
public:

    RigidBodyDemo(cWorld* world)
    {
        m_world = world;
        initialize();
    }

    ~RigidBodyDemo()
    {

    }

    void initialize();
    void step(double dt);
    void updateGraphics();
    void computeConstraints(double dt);

    cWorld* m_world;
    vector<RigidObject*> m_objects;
    Settings m_settings;
    CollisionDetector* m_collisionDetector;


};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H
