#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_PEGINHOLE_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_PEGINHOLE_H

#include "../../src/objects.hpp"

typedef enum {
    LOCAL_SOLVER,
    GLOBAL_SOLVER,
} SolverType;

struct Settings
{
    double K = 1000; double B = 100;
    double linearTolerance = 1e-6;
    int numLinearIterations = 20;
    double newtonTolerance = 1e-6;
    int numNewtonIterations = 20;
    double mass = 10;
    SolverType solverType = GLOBAL_SOLVER;
};

class PegInHole
{
public:

    PegInHole(cWorld* world, cGenericHapticDevicePtr device = NULL)
    {
        m_devicePtr = device; m_world = world;
        initialize();
    }

    ~PegInHole()
    {

    }

    void initialize(void);
    void step(double dt);
    void solve_local();
    void solve_global();

    void update_haptics();
    void update_graphics();

    cGenericHapticDevicePtr m_devicePtr; // the pointer to the haptic device
    cWorld* m_world;
    RigidObject* peg;
    RigidObject* block;
    Settings m_settings;

};
#endif //IMPLICITNONLINEARCOMPLEMENTARITY_PEGINHOLE_H
