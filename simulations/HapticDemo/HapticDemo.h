#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_HAPTICDEMO_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_HAPTICDEMO_H

#include "chai3d.h"
#include <Eigen/Dense>

using namespace chai3d;
using namespace std;
using namespace Eigen;

struct ColInfo
{
    Vector3d normal;
    Vector3d contactPt;
    double depth;
};

struct Plane
{
    cMesh* vis;
    Vector3d normal = Vector3d(0,0,1);
    Vector3d point = Vector3d(0,0,-0.1);

    Plane()
    {
        vis  = new cMesh();
        vis->m_material->setRed();
        cCreatePlane(vis,1,1,point);
    }
};

struct GodObject
{
    double s = 5.; // workspace scaling
    cShapeSphere* vis; // the visualizer of the haptic device
    double r = 0.025; // radius of haptic sphere
    cGenericHapticDevicePtr device; // the haptic device
    Vector3d x_d = Vector3d(0,0,0); // the position of the haptic device
    Vector3d xdot_d = Vector3d(0,0,0); // the velocity of the haptic device
    Vector3d x = Vector3d(0,0,0); // the god object position
    Vector3d x_tilde = Vector3d(0,0,0); //  the unconstrained object position
    Vector3d u = Vector3d(0,0,0); // the god object velocity
    Vector3d u_tilde = Vector3d(0,0,0); // the god object unconstrained velocity
    Matrix3d M = Matrix3d::Identity(); //  the god object mass matrix
    double K = 1000; // the coupling stiffness of the god object
    double B = 10; // the coupling damping of the god object

    GodObject()
    {
        vis  = new cShapeSphere(r);
    }

    Vector3d computeForce()
    {
        return K*(x - x_d);
    }

    void updateComplementarityProblem(double dt)
    {
        u_tilde = 0.5*(x_d - x) / dt ;
//        cout << "u_tilde = " << u_tilde.transpose() << endl;
        x_tilde = x_d;
        u.setZero();
    }

    void updateFromDevice()
    {
        cVector3d temp; device->getPosition(temp);
        x_d = s*temp.eigen();
        device->getLinearVelocity(temp);
        xdot_d = s*temp.eigen();
    }

    void updateVisualizer()
    {
        vis->setLocalPos(x);
    }

};

class HapticDemo
{
public:

    HapticDemo(cWorld* world, cGenericHapticDevicePtr a_device = NULL)
    {
        m_world = world;
        godObject = new GodObject;
        godObject->device = a_device;
        world->addChild(godObject->vis);
        plane = new Plane;
        world->addChild(plane->vis);
    }

    ~HapticDemo()
    {

    }

    void step(double dt);
    void updateGraphics();
    void updateHaptics(Vector3d& f,double dt);

    // the world members
    cWorld* m_world;
    GodObject* godObject;
    Plane* plane;

    // Simulation constants
    int numNewtonIt = 10;
    int backtrackingIt = 10;
    double alpha = 0.5; double beta = 0.2;
    double contactTol = 1e-6;

    // System submatrices
    MatrixXd C_h; // compliance matrix used for haptic feedback
    MatrixXd C_delassus; // the delassus operator
    MatrixXd J; // the constraint jacobian
    MatrixXd M; // the mass matrixc
    ColInfo* colInfo; // the collision information


};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_HAPTICDEMO_H
