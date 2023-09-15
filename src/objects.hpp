#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP

#include <Eigen/Dense>
#include <vector>
#include "chai3d.h"
#include "helper.hpp"
#include "constraints.hpp"

using namespace chai3d;
using namespace Eigen;
using namespace std;

class Contact;

enum  ObjectType {DEFORMABLE,RIGID};

class Object : public cMultiMesh
{
public:

    Object(int idx = -1, const std::string& meshFilename = "")
    {
        m_idx = idx;
        loadMesh(meshFilename);
    }
    ~Object() = default;

    virtual void timestep(double dt){}
    void loadMesh(const std::string& filename);
    virtual void set_local_pos(Vector3d pos){}
    virtual void set_local_rot(Quaterniond rot){}
    void initializeVisualizer();
    void scaleMesh(const double s);
    void set_is_static(bool is_static_){is_static = is_static_;}

    // the index
    int m_idx;  // auxiliary variable for global indexing
    double mass = 1.; // Mass.
    ObjectType type; // the type of objects
    vector<Contact*> m_contacts; // the contacts involving this object

    MatrixXi m_triangles; //  the triangles representing this object
    MatrixXd m_vertices; // the vertices of this object in the zero configuration
    MatrixXd m_normals; //  the normals for this object
    bool is_static = false; // is this object static
};

class RigidObject : public Object
{

public:

    RigidObject(int idx = -1, const std::string& meshFilename = "") : Object(idx,meshFilename)
    {
        type = RIGID;
    }

    MatrixXd M(void){return mass*Matrix3d::Identity();}
    VectorXd get_configuration(void);
    VectorXd get_configuration_unconstrained(void);
    virtual void set_local_pos(Vector3d pos) override;
    virtual void set_local_rot(Quaterniond rot) override;
    const MatrixXd& vertices(void){return m_vertices;}
    const MatrixXi& triangles(void){return m_triangles;}
    MatrixXd kinematic_map_G();
    static MatrixXd kinematic_map_G(Quaterniond);
    VectorXd generalized_pos();
    VectorXd generalized_vel();
    MatrixXd kinematic_map_H();
    MatrixXd generalized_mass();
    MatrixXd generalized_mass_inverse();
    void update_inertia_matrix();
    void compute_inertia_matrix();
    void update_mesh_position();

    Matrix3d I, Iinv = Matrix3d::Identity(); // Inertia and inverse inertia matrix (global)
    Matrix3d Ibody, IbodyInv = Matrix3d::Identity(); // Inertia and inverse inertia in the local body frame.
    Vector3d x = Vector3d::Zero(); // Position.
    Vector3d x_tilde = Vector3d::Zero(); // Unconstrained Position
    Quaterniond q = Quaterniond::Identity(); // Rotation
    Quaterniond q_tilde = Quaterniond::Identity(); // Unconstrained rotation.

    Matrix3d R = Matrix3d::Identity(); // Rotation matrix (auxiliary).
    Vector3d xdot = Vector3d::Zero(); // Linear velocity.
    Vector3d xdot_tilde = Vector3d::Zero(); // Unconstrained linear velocity
    Vector3d omega = Vector3d::Zero(); // Angular velocity.
    Vector3d omega_tilde = Vector3d::Zero(); // Unconstrained angular velocity.
    Vector3d f = Vector3d::Zero(); // Linear force.
    Vector3d tau = Vector3d::Zero(); // Angular force (torque).

    Vector3d fc = Vector3d::Zero(); // Linear constraint force.
    Vector3d tauc = Vector3d::Zero(); // Angular constraint force
};

class DeformableObject : public Object
{
public:

    DeformableObject(int idx = -1, const std::string& meshFilename = "") : Object(idx,meshFilename)
    {
        char* filenamePtr = const_cast<char*>(meshFilename.c_str());
        create_tetrahedral_mesh(filenamePtr);
        type = DEFORMABLE;
    }

    int ndof(void){return m_numVerts;}
    const MatrixXd& vertices(void){return m_vertices;}
    const MatrixXi& tetrahedra(void){return m_tetrahedra;}
    const MatrixXi& triangles(void){return m_triangles;}
    bool create_tetrahedral_mesh(char* filename);
    void update_mesh_position(void);
    MatrixXd mass_matrix();
    MatrixXd inverse_mass_matrix();

    int m_numVerts;             // the number of vertices
    int m_numTets;              // the number of tetrahedra

    VectorXd x;                  // Node positions.
    VectorXd x_last;             // Last positions.

    VectorXd xdot;               // Linear velocities.
    VectorXd xdot_last;          // Last velocities.
    VectorXd f;                  // External force.
    VectorXd fc;                 // Linear constraint force.
    MatrixXi m_tetrahedra;       // The tetrahedra

    vector<FixedConstraint*> fixed_constraints; // the fixed point constraints
    vector<LinearHookeanConstraint*> linear_hookean_constraints; //  the linear constraints
    vector<LinearCorotationalConstraint*> corotational_constraints; // the corotational constraints

};

// these are the haptic objects
class GodObject : public RigidObject {

public:

    GodObject(cGenericHapticDevicePtr a_device = NULL, int idx = -1, const std::string &meshFilename = "", bool visualizeGodObject = false) : RigidObject(idx, meshFilename) {
        device = a_device;
        if (visualizeGodObject)
        {

        }
    }

    void updateComplementarityProblem(double dt)
    {
        xdot_tilde = (x_d - x);
        x_tilde = x_d ;
        xdot.setZero();
    }

    void updateFromDevice()
    {
        cVector3d temp; device->getPosition(temp);
        x_d = s*temp.eigen();
    }

    void loadGodObjectMesh()
    {

    }

    Vector3d x_d = Vector3d::Identity(); // The position of the haptic device
    Quaterniond q_d = Quaterniond::Identity(); // The orientation of the haptic device
    cMultiMesh* vis;
    double s = 15; // workspace scaling
    double r = 0.025; // radius of haptic sphere
    cGenericHapticDevicePtr device; // the haptic device
    double K = 1000; // the coupling stiffness of the god object
    double B = 10; // the coupling damping of the god object
    cMultiMesh* godObjectVis; // visualize the god object
    double t = 0;

};

class DeformableGodObject : public DeformableObject
{
public:

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
