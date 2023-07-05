#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP

#include <Eigen/Dense>
#include <vector>
#include "chai3d.h"
#include "helper.hpp"

using namespace chai3d;
using namespace Eigen;
using namespace std;

class Contact;

enum  ObjectType {DEFORMABLE,RIGID};

class Object : public cMultiMesh
{
public:

    Object(int idx = -1)
    {
        m_idx = idx;
    }
    ~Object() = default;

    virtual void timestep(double dt){}
    virtual void set_local_pos(Vector3d pos){}
    virtual void set_local_rot(Quaterniond rot){}
    void scale(double s)
    {
        m_vertices *= s; this->getMesh(0)->scale(s);
    }

    // the index
    int m_idx;  // auxiliary variable for global indexing
    double mass = 1.; // Mass.
    ObjectType type; // the type of objects
    vector<Contact*> m_contacts; // the contacts involving this object

    //
    MatrixXi m_triangles; //  the triangles representing this object
    MatrixXd m_vertices; // the vertices of this object in the zero configuration
    MatrixXd m_normals; //  the normals for this object

};

class RigidObject : public Object
{

public:

    RigidObject(int idx = -1, const std::string& meshFilename = "") : Object(idx)
    {
        if (!this->loadFromFile(meshFilename))
        {
            cout << "COULD NOT LOAD MESH FILE" << endl;
        }
        x.setZero();
        x_last.setZero();
        q = Quaterniond::Identity();
        q_last = Quaterniond ::Identity();
        xdot.setZero();
        xdot_last.setZero();
        omega.setZero();
        omega_last.setZero();
        Ibody.setIdentity();
        IbodyInv.setIdentity();
        type = RIGID;
    }

    virtual void set_local_pos(Vector3d pos) override;
    virtual void set_local_rot(Quaterniond rot) override;
    const MatrixXd& vertices(void){return m_vertices;}
    const MatrixXi& triangles(void){return m_triangles;}
    void set_is_static(bool is_static_){is_static = is_static_;}
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
    void import_mesh_data();

    Eigen::Matrix3d I, Iinv;            // Inertia and inverse inertia matrix (global)
    Eigen::Matrix3d Ibody, IbodyInv;    // Inertia and inverse inertia in the local body frame.
    Eigen::Vector3d x;                  // Position.
    Eigen::Vector3d x_last;             // Last position.
    Eigen::Quaterniond q;               // Rotation
    Eigen::Quaterniond q_last;          // Last rotation.

    Eigen::Matrix3d R;                  // Rotation matrix (auxiliary).
    Eigen::Vector3d xdot;               // Linear velocity.
    Eigen::Vector3d xdot_last;          // Last velocity.
    Eigen::Vector3d omega;              // Angular velocity.
    Eigen::Vector3d omega_last;         // Last angular velocity.
    Eigen::Vector3d f;                  // Linear force.
    Eigen::Vector3d tau;                // Angular force (torque).

    Eigen::Vector3d fc;                 // Linear constraint force.
    Eigen::Vector3d tauc;               // Angular constraint force

    bool is_static = false;
};

class DeformableObject : public Object
{
public:

    DeformableObject(int idx = -1, const std::string& meshFilename = "") : Object(idx)
    {
        set_material_stiffness_matrix(9e9,0.5);
        char* filenamePtr = const_cast<char*>(meshFilename.c_str());
        create_tetrahedral_mesh(filenamePtr);
        type = DEFORMABLE;
    }

    int ndof(void){return m_numVerts;}
    const MatrixXd& vertices(void){return m_vertices;}
    const MatrixXi& tetrahedra(void){return m_tetrahedra;}
    const MatrixXi& triangles(void){return m_triangles;}
    void timestep(double dt)
    {

    }
    void set_material_stiffness_matrix(double,double);
    bool create_tetrahedral_mesh(char* filename);
    void compute_rest_volumes(void);
    void compute_rest_deformation_gradient(void);
    void compute_elasticity_matrix(void);
    void update_mesh_position(void);
    MatrixXd mass_matrix();
    MatrixXd inverse_mass_matrix();

    //
    int m_numVerts;
    int m_numTets;

    //
    MatrixXd K; // the material stiffness matrix
    MatrixXd E; // the elasticity matrix
    VectorXd volume_0;
    MatrixXd gradient_0;

    //
    MatrixXi m_tetrahedra;

};

class GodObject : public RigidObject
{
public:
};

class DeformableGodObject : public DeformableObject
{
public:

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
