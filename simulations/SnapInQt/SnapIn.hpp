
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_SNAPIN_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_SNAPIN_H

#include "../../src/objects.hpp"

typedef enum {
    LOCAL_SOLVER,
    GLOBAL_SOLVER,
} SolverType;

struct Settings
{
    // Haptic Settings
    double scaleFactor = 1.;
    double K = 1000; double B = 100;
    // Simulation settings
    double linearTolerance = 1e-6;
    int numLinearIterations = 20;
    double newtonTolerance = 1e-6;
    int numNewtonIterations = 20;
    double E = 1e10;
    double nu = 0.5;
    SolverType solverType = GLOBAL_SOLVER;
};

class Clip : public DeformableObject
{
public:
    Clip(){}
    ~Clip(){}

};

class Peg : public DeformableObject
{
public:
    Peg(){}
    ~Peg(){}
};

class SnapIn
{
public:

    SnapIn(cGenericHapticDevicePtr device = NULL)
    {
        m_devicePtr = device;
        initialize();
    }

    ~SnapIn()
    {

    }

    void initialize();
    void step(double dt);
    void solve_local();
    void solve_global();

    void update_haptics();
    void update_graphics();

    cGenericHapticDevicePtr m_devicePtr; // the pointer to the haptic device
    Clip* clip;
    Peg* peg;

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_SNAPIN_H
