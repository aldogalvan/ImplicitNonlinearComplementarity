#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP

#include <Eigen/Dense>
#include <vector>
#include "chai3d.h"
#include "helper.hpp"

using namespace chai3d;
using namespace Eigen;
using namespace std;

class Object : public cMultiMesh
{
public:

    Object(int idx, const std::string& meshFilename = "")
    {
        m_idx = idx;
        this->loadFromFile(meshFilename);
    }
    ~Object() = default;

    // the index
    int m_idx;

};

class RigidObject : public Object
{

public:

    RigidObject(int idx, const std::string& meshFilename = "") : Object(idx,meshFilename)
    {
        x.setZero();
        x_unconstrained.setZero();
        q = Quaterniond::Identity();
        q_unconstrained = Quaterniond ::Identity();
        xdot.setZero();
        omega.setZero();
        Ibody.setIdentity();
        IbodyInv.setIdentity();
    }

    void set_is_static(bool is_static_)
    {
        if (is_static_ == true){is_static = true;}
        else {is_static = false;}
    }
    void set_local_pos(Vector3d pos);
    MatrixXd kinematic_map_G();
    static MatrixXd kinematic_map_G(Quaterniond);
    MatrixXd kinematic_map_H();
    MatrixXd mass_matrix();
    MatrixXd inverse_mass_matrix();
    void update_inertia_matrix();
    void compute_inertia_matrix();
    void update_mesh_position();
    void import_mesh_data();

    double mass = 1.;                         // Mass.
    Eigen::Matrix3d I, Iinv;            // Inertia and inverse inertia matrix (global)
    Eigen::Matrix3d Ibody, IbodyInv;    // Inertia and inverse inertia in the local body frame.
    Eigen::Vector3d x;                  // Position.
    Eigen::Vector3d x_unconstrained;    // Position unconstrained
    Eigen::Quaterniond q;               // Orientation.
    Eigen::Quaterniond q_unconstrained; // Orientation unconstrained;

    Eigen::Vector3d xdot;               // Linear velocity.
    Eigen::Vector3d omega;              // Angular velocity.
    Eigen::Vector3d f;                  // Linear force.
    Eigen::Vector3d tau;                // Angular force (torque).

    Eigen::Vector3d fc;                 // Linear constraint force.
    Eigen::Vector3d tauc;               // Angular constraint force

    //
    MatrixXi triangles; //  the triangles representing this object
    MatrixXd vertices; // the vertices of this object in the zero configuration
    MatrixXd normals; //  the normals for this object

    bool is_static = false;
};

class GodObject : public RigidObject
{

};

class DeformableObject : public Object
{
public:

    DeformableObject(int idx) : Object(idx)
    {

    }

    void tetrahedralize();

    //
    MatrixXi triangles;
    MatrixXi tetrahedra;
    MatrixXd vertices;
    MatrixXd vertices_dot;
    Vector3d com;
};

class DeformableGodObject : public DeformableObject
{
public:


};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
