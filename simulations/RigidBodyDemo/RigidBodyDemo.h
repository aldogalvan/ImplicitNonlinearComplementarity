
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H

#include "chai3d.h"
#include <Eigen/Dense>
#include "../../src/objects.hpp"
#include "../../src/collision.hpp"
#include "../../src/helper.hpp"

struct Settings
{
    // Haptic Settings
    double scaleFactor = 1.;
    double K = 1000; double B = 100;
    // Simulation settings
    double linearTolerance = 1e-6;
    int numLinearIterations = 20;
    double newtonTolerance = 1e-6;
    int newtonIterations = 20;
    double E = 1e10;
    double nu = 0.5;
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
    VectorXd computeResidual(const VectorXd&, const Quaterniond&, const VectorXd &, const VectorXd&, const VectorXd&, const MatrixXd &,
                             const vector<Contact*>&, double);
    void backtrackingLineSearch(double &, double &, const VectorXd &, const VectorXd &, const VectorXd &,
                                const VectorXd &,const VectorXd&, const VectorXd&, const MatrixXd &,
                                const vector<Contact*>&, double, double, double, int);

    cWorld* m_world;
    vector<RigidObject*> m_objects;
    Settings m_settings;
    CollisionDetector* m_collisionDetector;
};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_RIGIDBODYDEMO_H
