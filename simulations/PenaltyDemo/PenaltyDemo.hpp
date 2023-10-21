

#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H

#include "objects.hpp"
#include "collision.hpp"
#include "contact.hpp"

using namespace std;

class PenaltyDemo
{
public:

    PenaltyDemo(cWorld* world, cGenericHapticDevicePtr device = NULL)
    {
        m_devicePtr = device; m_world = world;
        initialize();
        vis = new cMultiMesh();
        vis->loadFromFile("/home/aldo/ImplicitNonlinearComplementarity/resources/PegInHole/peg.obj");
        vis->m_material->setYellowLightGoldenrod();
        vis->scale(0.025);
        vis->scale(0.99);
        vis->setWireMode(true);
        cout << "center = " << -vis->getMesh(0)->getCenterOfMass() << endl;
        vis->getMesh(0)->offsetVertices(-vis->getMesh(0)->getCenterOfMass());
        world->addChild(vis);
    }

    ~PenaltyDemo()
    {

    }

    void initialize(void);
    void step(double dt);
    void updateGraphics();
    void updateHaptics(Vector3d& f);
    VectorXd computeResidual(const VectorXd&, const VectorXd &,const VectorXd &, const VectorXd &,
                             const VectorXd&, const VectorXd&, const MatrixXd &,
                             const vector<Contact*>&, double);
    void backtrackingLineSearch(double &, double &, const VectorXd &, const VectorXd &, const VectorXd &,
                                const VectorXd &,const VectorXd&, const VectorXd&, const MatrixXd &,
                                const vector<Contact*>&, double, double, double, int);

    cGenericHapticDevicePtr m_devicePtr; // the pointer to the haptic device
    cWorld* m_world;
    GodObject* peg;
    RigidObject* block;
    cMultiMesh* vis;


    // Simulation constants
    int numNewtonIt = 10;
    int backtrackingIt = 10;
    double alpha = 0.5; double beta = 0.2;

    // System submatrices
    MatrixXd C_h; // compliance matrix used for haptic feedback
    MatrixXd J; // the constraint jacobian
    MatrixXd C; // the compliance matrix
    MatrixXd M; // the mass matrix

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_PENALTYDEMO_H
