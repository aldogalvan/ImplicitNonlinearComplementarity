#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP

#include <Eigen/Dense>
#include <vector>
#include <set>
#include "chai3d.h"
#include "helper.hpp"
#include "constraints.hpp"

using namespace chai3d;
using namespace Eigen;
using namespace std;

class Contact;
struct Sensor;

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
    virtual void scaleMesh(const double s);
    void set_is_static(bool is_static_){is_static = is_static_;}

    // the index
    int m_idx;  // auxiliary variable for global indexing
    double mass = 0.01; // Mass.
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
    virtual void scaleMesh(const double s) override;
    const MatrixXd& vertices(void){return m_vertices;}
    const MatrixXi& tetrahedra(void){return m_tetrahedra;}
    const MatrixXi& triangles(void){return m_triangles;}
    bool create_tetrahedral_mesh(char* filename);
    void update_mesh_position(void);
    MatrixXd& mass_matrix();
    MatrixXd& inverse_mass_matrix();
    SparseMatrix<double>& sparse_mass_matrix();
    SparseMatrix<double>& sparse_inverse_mass_matrix();

    int m_numVerts;              // the number of vertices
    int m_numTets;               // the number of tetrahedra

    VectorXd x;                  // Node positions.
    VectorXd x_tilde;            // Unconstrained positions

    VectorXd xdot;               // Linear velocities.
    VectorXd xdot_tilde;         // Unconstrained linear velocities velocities.
    VectorXd f;                  // External force.
    VectorXd fc;                 // Linear constraint force.
    MatrixXi m_tetrahedra;       // The tetrahedra

    MatrixXd M;                  // the mass matrix
    MatrixXd M_inv;              // the inverse mass matrix
    SparseMatrix<double> M_sparse;      // the sparse mass matrix
    SparseMatrix<double> M_inv_sparse;  // the sparse inverse mass matrix

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
        xdot_tilde = 0.5*(x_d - x) / dt;
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


class PenaltyObject : public RigidObject
{

public:

    PenaltyObject( bool a_useVis = true, int idx = -1, const std::string& meshFilename = "") : RigidObject(idx,meshFilename)
    {
        visualizeSensors = a_useVis;
    }

    ~PenaltyObject()
    {

    }

    void wrapMesh();
    void updateSensors();

    Vector3d x_minus_1 = Vector3d::Zero(); // Position.
    Quaterniond q_minus_1 = Quaterniond::Identity(); // Rotation
    Vector3d xdot_minus_1 = Vector3d::Zero(); // Linear velocity.
    Vector3d omega_minus_1 = Vector3d::Zero(); // Angular velocity.

    vector<set<int>> vertexNeighbors;
    vector<Vector3d> vertexNormals;
    vector<vector<Sensor*>> m_sensors;
    bool visualizeSensors;
    // the visualization of the rays in this world
    vector<cShapeLine*> m_rayViz;

    double density = 1e5;
};

class PenaltyGodObject : public PenaltyObject
{
public:

    PenaltyGodObject(cGenericHapticDevicePtr a_device = NULL, bool a_useVis = true, int idx = -1, const std::string& meshFilename = "") : PenaltyObject(a_useVis,idx,meshFilename)
    {
        device = a_device;
    }

    ~PenaltyGodObject(){}

    void updateFromDevice()
    {
        cVector3d temp;
        device->getPosition(temp);
        cMatrix3d R;
        device->getRotation(R);
        q_d = Quaterniond(R.eigen());
        x_d = s*temp.eigen();
        device->getLinearVelocity(temp);
        xdot_d = temp.eigen();
        device->getAngularVelocity(temp);
        omega_d = temp.eigen();

    }

    Vector3d x_d = Vector3d::Identity(); // The position of the haptic device
    Quaterniond q_d = Quaterniond::Identity(); // The orientation of the haptic device
    Vector3d xdot_d = Vector3d::Zero();
    Vector3d omega_d = Vector3d::Zero();
    cMultiMesh* vis;
    bool m_recordTrajectory = false;
    bool m_readTrajectory = false;
    double s = 15; // workspace scaling
    double r = 0.025; // radius of haptic sphere
    cGenericHapticDevicePtr device; // the haptic device
    double Kc = 200; // the coupling stiffness of the god object
    double Bc = 1; // the coupling damping of the god object
    double Kt = 1000; // the coupling stiffness of the god object
    double Bt = 1; // the coupling damping of the god object
    cMultiMesh* godObjectVis; // visualize the god object
    double t = 0;
};

struct Sensor
{
    Vector3d globalStartPos; // the global start position
    Vector3d globalNormal; // the global normal
    Vector3d startPos; // start position defined in local frame of reference
    Vector3d normal; // end position defined in local frame of reference
    double len = 0.01;
    double k = 0; double b = 0;
    PenaltyObject* m_parent;
};


#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
