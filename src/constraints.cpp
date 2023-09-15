#include "constraints.hpp"

Matrix<double, 9, 12> computedFdx(const Eigen::Matrix3d& DmInv) {
    const double m = DmInv(0, 0);
    const double n = DmInv(0, 1);
    const double o = DmInv(0, 2);
    const double p = DmInv(1, 0);
    const double q = DmInv(1, 1);
    const double r = DmInv(1, 2);
    const double s = DmInv(2, 0);
    const double t = DmInv(2, 1);
    const double u = DmInv(2, 2);

    Eigen::Matrix<double, 9, 12> PFPu = Eigen::Matrix<double, 9, 12>::Zero();
    PFPu.row(0) << -m - p - s, 0.0, 0.0, m, 0.0, 0.0, p, 0.0, 0.0, s, 0.0, 0.0;
    PFPu.row(1) << 0.0, -n - q - t, 0.0, 0.0, n, 0.0, 0.0, q, 0.0, 0.0, t, 0.0;
    PFPu.row(2) << 0.0, 0.0, -o - r - u, 0.0, 0.0, o, 0.0, 0.0, r, 0.0, 0.0, u;
    PFPu.row(3) << -n - q - t, n, 0.0, q, 0.0, 0.0, t, 0.0, 0.0, 0.0, 0.0, 0.0;
    PFPu.row(4) << m, 0.0, 0.0, -m, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    PFPu.row(5) << 0.0, -n - q - t, 0.0, 0.0, n, 0.0, 0.0, q, 0.0, 0.0, t, 0.0;
    PFPu.row(6) << -o - r - u, 0.0, o, r, 0.0, 0.0, u, 0.0, 0.0, 0.0, 0.0, 0.0;
    PFPu.row(7) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
    PFPu.row(8) << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;

    return PFPu;
}


//////////////////////////////////////////////////////////////////////
////////////////////// THE FIXED CONSTRAINT //////////////////////////
//////////////////////////////////////////////////////////////////////

VectorXd FixedConstraint::computeJacobian(VectorXd &x)
{
    J.resize(3,1);
    auto p = x.block<3,1>(3*v0,0);
    J = (p0 - p).normalized();
    return J;
}

double FixedConstraint::compute_Cb(VectorXd &x)
{
    auto p = x.block<3,1>(3*v0,0);
    Cb = (p0 - p).norm();
    return Cb;
}

//////////////////////////////////////////////////////////////////////
///////////////////// THE LINEAR CONSTRAINT //////////////////////////
//////////////////////////////////////////////////////////////////////

void LinearHookeanConstraint::computeC()
{
    C.resize(6,6);
    auto K_inv = K.inverse();
    C = K_inv / V0;
}

void LinearHookeanConstraint::compute_Dm(VectorXd& x)
{
    MatrixXd Dm(3,3); Dm.setZero();
    Dm_inv.resize(3,3); Dm_inv.setZero();

    Vector3d p0 = x.block<3,1>(v0*3,0);
    Vector3d p1 = x.block<3,1>(v1*3 ,0);
    Vector3d p2 = x.block<3,1>(v2*3 ,0);
    Vector3d p3 = x.block<3,1>(v3*3,0);

    Vector3d o1 = p1 - p0; Vector3d o2 = p2 - p0; Vector3d o3 = p3 - p0;
    Dm.col(0) = o1; Dm.col(1) = o2; Dm.col(2) = o3;
    Dm_inv= Dm.inverse();
}


MatrixXd LinearHookeanConstraint::compute_F(VectorXd& x)
{
    MatrixXd F = Matrix3d::Zero();
    Matrix3d Ds;

    Vector3d p0 = x.block<3,1>(v0*3,0);
    Vector3d p1 = x.block<3,1>(v1*3 ,0);
    Vector3d p2 = x.block<3,1>(v2*3 ,0);
    Vector3d p3 = x.block<3,1>(v3*3,0);

    Vector3d o1 = p1 - p0; Vector3d o2 = p2 - p0; Vector3d o3 = p3 - p0;
    Ds.col(0) = o1; Ds.col(1) = o2; Ds.col(2) = o3;
    F = Ds*Dm_inv;
    return F;
}

void LinearHookeanConstraint::computeRestVolume(VectorXd& x)
{
    Vector3d p0 = x.block<3,1>(v0*3,0);
    Vector3d p1 = x.block<3,1>(v1*3,0);
    Vector3d p2 = x.block<3,1>(v2*3,0);
    Vector3d p3 = x.block<3,1>(v3*3,0);

    V0 = computeTetrahedronVolume(p0,p1,p2,p3);
}

void LinearHookeanConstraint::setMaterialStiffnessMatrix(double E, double nu)
{
    // Compute the Lame parameters
    double lambda = (E * nu) / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    // Create the 6x6 material stiffness matrix (assuming 3D)
    Eigen::MatrixXd stiffnessMatrix(6, 6);

    stiffnessMatrix << lambda + 2 * mu, lambda, lambda, 0, 0, 0,
            lambda, lambda + 2 * mu, lambda, 0, 0, 0,
            lambda, lambda, lambda + 2 * mu, 0, 0, 0,
            0, 0, 0, mu, 0, 0,
            0, 0, 0, 0, mu, 0,
            0, 0, 0, 0, 0, mu;
    K = stiffnessMatrix;
}

MatrixXd LinearHookeanConstraint::computeGreenStrain(MatrixXd& F)
{
    return 0.5*(F.transpose()*F - Matrix3d::Identity());
}

MatrixXd LinearHookeanConstraint::computeJacobian(VectorXd &x)
{
    MatrixXd F = compute_F(x);
    Matrix3d dSdF = 0.5 * (F + F.transpose());
    MatrixXd dFdx = computedFdx(Dm_inv);
    MatrixXd dSdx(6,12);

    for (int i = 0 ; i < 12 ; i++)
    {
        dSdx.col(i) = strainTensorToVector(dSdF * reshape(dFdx.col(i)));
    }

    return dSdx;
}

VectorXd LinearHookeanConstraint::compute_Cb(Eigen::VectorXd &x)
{
    MatrixXd F = compute_F(x);
    VectorXd s = strainTensorToVector(0.5 * (F*F.transpose() - Matrix3d::Identity()));
    return s;
}


//////////////////////////////////////////////////////////////////////
/////////////////// THE COROTATIONAL CONSTRAINT //////////////////////
//////////////////////////////////////////////////////////////////////
