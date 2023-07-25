#include <Eigen/Dense>
#include "helper.hpp"

using namespace Eigen;

class FixedConstraint
{
public:

    FixedConstraint(int v0_ , Vector3d p0_)
    {
        v0 = v0_;
        p0 = p0_;
    }

    void updateConstraint(MatrixXd& x);
    void computeJacobian(MatrixXd& x); // computes the jacobian wrt to vertex
    void compute_Cb(MatrixXd& x);

    int v0; // the vertex this constraint acts on
    Vector3d p0; // the fixed point position
    MatrixXd J; // the jacobian
    double Cb; // the constraint value

};

class LinearHookeanConstraint
{
public:

    LinearHookeanConstraint(int v0_, int v1_, int v2_, int v3_, int tidx_, VectorXd& x,
                     double E = 9e9, double nu = 0.48)
    {
        v0 = v0_; v1 = v1_; v2 = v2_; v3 = v3_;
        tidx = tidx_;

        setMaterialStiffnessMatrix(E,nu);
        computeRestVolume(x);
        computeC(x);
        compute_Dm(x);
    }

    void updateConstraint();
    void computeC(VectorXd& x);
    void compute_Dm(VectorXd& x);
    MatrixXd compute_F(VectorXd& x);
    void computeRestVolume(VectorXd& x);
    void setMaterialStiffnessMatrix(double E, double nu);
    MatrixXd computeGreenStrain(MatrixXd& F);
    VectorXd compute_Cb(VectorXd& x);
    MatrixXd computeJacobian(VectorXd& x);
    MatrixXd PK1(MatrixXd& F);

public:

    MatrixXd K; // the material stiffness matrix (9x9)
    MatrixXd E; // the strain tensor
    MatrixXd C; // the compliance matrix (9x9)
    MatrixXd Dm_inv;
    double V0; //  the rest volumes
    int tidx; // the tetrahedron idx
    int v0, v1, v2, v3; // the vertex indices


};

class LinearCorotationalConstraint
{
public:

    LinearCorotationalConstraint(int v0_, int v1_, int v2_, int v3_, int tidx_, MatrixXd& x)
    {
        v0 = v0_; v1 = v1_; v2 = v2_; v3 = v3_;
        tidx = tidx_;
    }

    int tidx;
    int v0, v1, v2, v3; // the vertex indices
    MatrixXd K; // the material stiffness matrix
    MatrixXd E; // the strain tensor
    MatrixXd C; // the compliance matrix
    MatrixXd Dm;
    MatrixXd Dm_inv;
    MatrixXd F; // the deformation gradient
    double V0; //  the rest volumes
    MatrixXd R; // the initial rotation
};