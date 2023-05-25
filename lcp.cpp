#include "lcp.hpp"
#include <drake/solvers/moby_lcp_solver.h>

// Function to solve the LCP using the Lemke algorithm
inline bool lemkeSolver(const Eigen::MatrixXd& M, const Eigen::VectorXd& q, Eigen::VectorXd& z) {
    const int n = q.size();
    const int m = M.rows();

    // Set up the extended matrix for Lemke's algorithm
    Eigen::MatrixXd L(n + m, n + m);
    Eigen::VectorXd b(n + m);
    L.block(0, 0, n, n) = M;
    L.block(0, n, n, n).setIdentity();
    L.block(n, 0, m, n) = -M.transpose();
    L.block(n, n, m, n).setZero();
    b.segment(0, n) = -q;
    b.segment(n, m).setZero();

    // Initialize the basic variables
    Eigen::VectorXi basic(n + m);
    for (int i = 0; i < n + m; ++i) {
        basic(i) = i;
    }

    // Lemke's algorithm
    int leaving = -1;
    while (true) {
        // Check if the current solution is valid
        bool valid = true;
        for (int i = 0; i < n; ++i) {
            if (z(i) < 0.0) {
                valid = false;
                break;
            }
        }
        if (valid) {
            return true;  // Valid solution found
        }

        // Determine the entering variable
        int entering = -1;
        for (int i = 0; i < n + m; ++i) {
            if (z(i) < 0.0) {
                entering = i;
                break;
            }
        }

        if (entering == -1) {
            return false;  // No entering variable found
        }

        // Solve for the leaving variable
        leaving = -1;
        for (int i = 0; i < n + m; ++i) {
            if (L(i, entering) > 0.0) {
                if (leaving == -1 || z(i) / L(i, entering) < z(leaving) / L(leaving, entering)) {
                    leaving = i;
                }
            }
        }
        if (leaving == -1) {
            return false;  // Unbounded LCP
        }

        // Perform the pivot operation
        double pivot = L(leaving, entering);
        for (int i = 0; i < n + m; ++i) {
            if (i == leaving) {
                continue;
            }
            double ratio = L(i, entering) / pivot;
            for (int j = 0; j < n + m; ++j) {
                L(i, j) -= ratio * L(leaving, j);
            }
            b(i) -= ratio * b(leaving);
        }

        // Update the basic variables
        basic(leaving) = entering;

        // Update the solution
        z(leaving) = b(leaving) / pivot;

        // Update z for non-basic variables
        for (int i = 0; i < n + m; ++i) {
            if (i != leaving) {
                z(i) -= L(i, entering) * z(leaving) / pivot;
            }
        }
    }
}

Eigen::MatrixXd getDirectionVectors(const Eigen::Vector3d& N) {
    Eigen::MatrixXd D(3, 4);

    // Normalize the normal vector
    Eigen::Vector3d normalizedN = N.normalized();

    // Compute the first direction vector using cross product
    Eigen::Vector3d V1;
    V1(0) = -normalizedN(1);
    V1(1) = normalizedN(0);
    V1(2) = 0;
    V1.normalize();

    auto dot = V1.dot(normalizedN);
    assert(abs(dot) < 1e-6);

    // Compute the second direction vector using cross product
    Eigen::Vector3d V2 = normalizedN.cross(V1).normalized();

    // Store the direction vectors in the matrix D
    D.col(0) = V1;
    D.col(1) = -V1;
    D.col(2) = V2;
    D.col(3) = -V2;

    return D;
}

MatrixXd contactJacobian(Vector3d N,Vector3d R)
{
    // three dof
    MatrixXd J(3,3);
    J.setIdentity();
    // TODO: Six dof

    return J;
}

bool LCP::setup_lcp(VectorXd& tau, vector<ColInfo*>& collisions)
{

    cout << collisions.size()<< endl;
    if (collisions.size() > 10)
        return 0;

    int numCollisions = collisions.size();
    // assume its only two objects colliding
    int numBodies = 2;
    int numDOF = 3*numBodies;
    double dt = 0.001;
    int basis_dim = 4;

    assert(tau.size() == numDOF);

    // generalized mass matrix
    M_inv.resize(6,6); M_inv.setIdentity();

    MatrixXd N(numDOF,numCollisions);
    MatrixXd B(numDOF,basis_dim*numCollisions);// four basis vectors for friction
    MatrixXd E(numCollisions*basis_dim,numCollisions);
    MatrixXd Mu(numCollisions,numCollisions);
    Mu.setIdentity(); Mu *= mu_;

    // friction cone basis matrix
    for (int colIdx = 0; colIdx < numCollisions; colIdx ++)
    {
        Vector3d normal = collisions[colIdx]->normal;
        MatrixXd D = getDirectionVectors(normal);

        //
        int bodyIdxA = collisions[colIdx]->bodyIdxA;
        int bodyIdxB = collisions[colIdx]->bodyIdxB;

        // contact jacobians
        MatrixXd Ja = -contactJacobian(normal,Vector3d(0,0,0));
        MatrixXd Jb = contactJacobian(normal,Vector3d(0,0,0));

        // add to normal jacobian
        N.block<3,1>(3*bodyIdxA,colIdx) = Ja*normal;
        N.block<3,1>(3*bodyIdxB,colIdx) = Jb*normal;

        // add to tangential jacobian
        B.block<3,4>(3*bodyIdxA,basis_dim*colIdx) = Ja*D;
        B.block<3,4>(3*bodyIdxB,basis_dim*colIdx) = Jb*D;

        // fill e matrix
        E.block<4,1>(basis_dim*colIdx,colIdx) = Vector4d(1,1,1,1);
    }

    // A Matrix
    A.resize(numCollisions + numCollisions*basis_dim + numCollisions,numCollisions + numCollisions*basis_dim + numCollisions);
    A.block(0,0,numCollisions,numCollisions) = dt*N.transpose()*M_inv*N;
    A.block(0,numCollisions,numCollisions,numCollisions*basis_dim) = dt*N.transpose()*M_inv*B;
    A.block(0,numCollisions*basis_dim + numCollisions,numCollisions,numCollisions) = MatrixXd(numCollisions,numCollisions).setZero();
    A.block(numCollisions,0,basis_dim*numCollisions,numCollisions) = dt*B.transpose()*M_inv*N;
    A.block(0,numCollisions,basis_dim*numCollisions,numCollisions*basis_dim) = dt*B.transpose()*M_inv*B;
    A.block(0,numCollisions,basis_dim*numCollisions,numCollisions) = E;
    A.block(numCollisions + basis_dim*numCollisions,0,numCollisions,numCollisions) = Mu;
    A.block(numCollisions + basis_dim*numCollisions,numCollisions,numCollisions,numCollisions*basis_dim) = E.transpose();
    A.block(numCollisions + basis_dim*numCollisions,numCollisions*basis_dim + numCollisions,numCollisions,numCollisions) = MatrixXd(numCollisions,numCollisions).setZero();

    cout << "New Problem" << endl << endl;
    cout << A << endl << endl;
    cout << b.transpose() << endl << endl;
    // b vector
    b.resize(numCollisions + numCollisions*basis_dim + numCollisions);
    b.head(numCollisions) = N.transpose()*M_inv*tau;
    b.segment(numCollisions,numCollisions*basis_dim ) = B.transpose()*M_inv*tau;
    b.tail(numCollisions) = VectorXd(numCollisions).setZero();

    return 1;
}

void LCP::solve_lcp_lemke(Eigen::VectorXd *lambda)
{
    //lemkeSolver(A,b,lambda);

    // Create an instance of the MobyLCPSolver class
    drake::solvers::MobyLCPSolver<double> lcp_solver;
    lcp_solver.SolveLcpFast(A,b,lambda);
}