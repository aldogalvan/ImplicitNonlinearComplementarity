
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_BEAMDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_BEAMDEMO_H

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
    int newtonIterations = 10;
    double E = 1e10;
    double nu = 0.5;
};

class BeamDemo
{
public:

    BeamDemo(cWorld* world)
    {
        m_world = world;
        initialize();
    }

    ~BeamDemo()
    {

    }

    void initialize();
    void step(double dt);
    void updateGraphics();
    void computeConstraints(double dt);

    cWorld* m_world;
    vector<DeformableObject*> m_objects;
    Settings m_settings;
    CollisionDetector* m_collisionDetector;

    // Simulation matrices
    MatrixXd M; // mass matrix
    MatrixXd M_inv; // inverse mass matrix
    VectorXd u; // the velocity vector
    VectorXd u_tilde; // the unconstrained velocity vector


};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_BEAMDEMO_H
