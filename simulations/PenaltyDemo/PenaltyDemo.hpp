

#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H

#include "objects.hpp"
#include "collision.hpp"
#include "contact.hpp"
#include <fstream>
#include <Eigen/Dense>

using namespace std;

struct CollisionResult {
    bool collision;
    double distance; // the penetration distance
    Eigen::Vector3d p0; // absolute contact point
    Vector3d normal; // the normal of the contact
    Sensor* sensor; // pointer to the sensor

};

class PenaltyDemo
{
public:

    PenaltyDemo(cWorld* world, cGenericHapticDevicePtr device = NULL)
    {
        m_devicePtr = device; m_world = world;
        initialize();
    }

    ~PenaltyDemo()
    {

    }

    void initialize(void);
    void step(double dt);    void updateGraphics();
    void updateHaptics(Vector3d& f);
    void linearizedBDF1Solver(double, vector<Collision>);
    void BDF1Solver(double h, vector<Collision>);
    void linearizedBDF2Solver(double h, vector<Collision>);
    void BDF2Solver(double h,vector<Collision>);

    vector<CollisionResult> computeCollisions(vector<Collision> potentialCollisions);
    cGenericHapticDevicePtr m_devicePtr; // the pointer to the haptic device
    cWorld* m_world;
    PenaltyGodObject* spoon;
    PenaltyObject* mug;

    // Simulation constants
    double simTime = 0;
    int maxIt = 20;
    bool useContact = 1;

    // System submatrices
    MatrixXd C_h; // compliance matrix used for haptic feedback
    MatrixXd J; // the constraint jacobian
    MatrixXd C; // the compliance matrix
    MatrixXd M; // the mass matrix
    bool m_recordTrajectory = false;
    bool m_readTrajectory = false;
    bool m_recordData = false;
    double h_minus_1 = 0.001;
    std::ifstream* inputFile;
    std::ofstream* outputFile;
    std::ofstream* dataFile;


};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H
