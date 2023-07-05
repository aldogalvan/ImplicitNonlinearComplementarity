#include "lcp.hpp"
#include <drake/solvers/moby_lcp_solver.h>
#include <algorithm>

class LineSearchMonitor {

private:
    double alpha;
    double beta;
    double tau_min;
    int line_search_iter;
    double tau;
    double fk;
    double f0;
    double grad_f0;

public:
    LineSearchMonitor(double alpha = 0.5, double beta = 0.0001, double tau_min = std::numeric_limits<double>::epsilon())
            : alpha(alpha), beta(beta), tau_min(tau_min), line_search_iter(0), tau(1.0),
              fk(std::numeric_limits<double>::max()), f0(std::numeric_limits<double>::max()), grad_f0(std::numeric_limits<double>::max()) {}

    void init(double f0, double grad_f0){
        this->f0 = f0;
        this->grad_f0 = this->beta * grad_f0;
    }

    void reset() {
        line_search_iter = 0;
        tau = 1.0;
    }

    bool finished(double fk) {
        this->fk = fk;
        return step_length_found() || check_step_length();
    }

    bool step_length_found() {
        if (fk <= f0 + tau * grad_f0) {
            return true;
        } else {
            tau *= alpha;
            return false;
        }
    }

    bool check_step_length() {
        return tau * tau < tau_min;
    }

    void step() {
        line_search_iter++;
    }

    double step_length() {
        return tau;
    }

    int line_search_steps() {
        return line_search_iter;
    }
};

// Function to solve the LCP using the Lemke algorithm
inline bool lemkeSolver(const MatrixXd& M, const VectorXd& q, VectorXd& z) {

    const int n = q.size();
    const int m = M.rows();

    // Check if z has the appropriate size
    if (z.size() != n + m) {
        z.resize(n + m);
        z.setZero();
    }

    // Set up the extended matrix for Lemke's algorithm
    MatrixXd L(n + m, n + m);
    VectorXd b(n + m);
    L.block(0, 0, n, n) = M;
    L.block(0, n, n, n).setIdentity();
    L.block(n, 0, m, n) = -M.transpose();
    L.block(n, n, m, n).setZero();
    b.segment(0, n) = -q;
    b.segment(n, m).setZero();

    // Initialize the basic variables
    VectorXi basic(n + m);
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

double compute_objective(const VectorXd& phi)
{
    return 0.5*phi.transpose()*phi;
}

VectorXd fischer(const VectorXd& y, const VectorXd& x)
{
    return sqrt(x.array().square() + y.array().square()) - x.array() - y.array();
}



void fischer_newton_solver(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, int max_iter = 0, double tol_rel = 0.00001,
                           double tol_abs = std::numeric_limits<double>::epsilon() * 10, string solver = "perturbation",
                           bool profile = true, LineSearchMonitor searchstrategymonitor = LineSearchMonitor())
{
    std::map<int, std::string> msg = {
            {1, "preprocessing"},
            {2, "iterating"},
            {3, "relative"},
            {4, "absolute"},
            {5, "stagnation"},
            {6, "local minima"},
            {7, "nondescent"},
            {8, "maxlimit"},
    };

    int N = b.size();
    int flag = 1;

    assert(x.size() == N && x.cols() == 1);
    assert(A.rows() == N && A.cols() == N);
    assert(b.size() == N && b.cols() == 1);

    if (max_iter == 0)
        max_iter = std::floor(N / 2.0);

    max_iter = std::max(max_iter, 1);

    const double h = 1e-7;
    const double gamma = 1e-28;
    const double eps = std::numeric_limits<double>::epsilon();
    const double rho = std::numeric_limits<double>::epsilon();
    const double gmres_tol = 10 * eps;

    bool warm_start = false;
    int max_warm_iterate = 5;
    int max_bgs_iterate = 5;

//    if (warm_start)
//    {
//        x = fischer_warm_start(A, x, b, max_warm_iterate, max_bgs_iterate);
//    }

    Eigen::VectorXd convergence(max_iter + 1);
    double err = 1e20;
    int iterate = 1;
    flag = 2;
    int grad_steps = 0;
    int grad_steps_max = max_iter;
    bool use_grad_steps = true;

    while (iterate <= max_iter)
    {
        cout << "Iterate: " << iterate << endl;
        bool take_grad_step = false;

        auto y = A*x + b;
        assert(y.cols() == 1);

        // calculate the fischer function value
        auto phi = fischer(y, x);
        assert(y.cols() == 1);
        auto old_err  = err;
        // calculate merit value
        err = compute_objective(phi);

        if (profile)
            convergence[iterate - 1] = err;

        if (grad_steps > grad_steps_max)
            use_grad_steps = false;

        // test the stopping criterias used
        auto rel_err = abs(err - old_err)/abs(old_err);

        if (rel_err < tol_rel)
        {
            if (use_grad_steps) {
                take_grad_step == true;
            } else if (rel_err < eps) {
                flag = 3;
                break;
            }
        }

        if (err < tol_abs)
        {
            flag = 4;
            break;
        }

        // Solving the Newton System
        auto dx = VectorXd(N).setZero();

        // Create a temporary VectorXd object to store the intermediate results
        Eigen::VectorXd S = Eigen::VectorXd::Zero(phi.size());

        cout << "phi : " << phi.transpose() << endl;

        // Perform element-wise comparison and logical operations
        for (int i = 0; i < phi.size(); ++i) {
            if (std::abs(phi(i)) < gamma && std::abs(x(i)) < gamma) {
                S(i) = 1.0;
            }
        }

        // Create a temporary VectorXd object to store the intermediate results
        Eigen::VectorXd I = Eigen::VectorXd::Ones(S.size()) - S;

        // Find the indices where I is true (non-zero)
        std::vector<int> idx;
        for (int i = 0; i < I.size(); ++i) {
            if (I(i) != 0.0) {
                idx.push_back(i);
            }
        }

        auto px = VectorXd(x.size()).setZero();
        px = x;

        cout << "x : " << x.transpose() << endl;

        Eigen::VectorXd direction = x.array().sign();

        // Replace zeros in direction with ones
        direction = direction.array().select(1.0, direction);

        //! START PERTURBATION
        // Update px based on S, gamma, and direction
        cout << "gamma : " << gamma << endl;
        cout << "direction : " << direction.transpose() << endl;
        cout << "S : " << S.transpose() << endl;
        px = gamma * direction.array() * S.array().cast<double>();
        cout << "px : " << px.transpose() << endl;
        cout << "y : " << y.transpose() << endl;
        auto p = VectorXd(px.size()).setZero();
        p = px.array() / (sqrt(y.array().square() + px.array().square())) - 1;
        cout << "p : " << p.transpose() << endl;
        assert(p.size() == N);
        auto q = VectorXd(y.size()).setZero();
        q = y.array() / (sqrt(y.array().square() + px.array().square())) - 1;
        cout << "q : " << q.transpose() << endl;
        assert(q.size() == N);
        MatrixXd J = p.asDiagonal() * MatrixXd::Identity(N,N) + q.asDiagonal()*A;
        assert(J.cols() == N && J.rows() == N);
        // Compute the condition number of J
        Eigen::JacobiSVD<Eigen::MatrixXd> svd(J);
        double cond = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);

        // Compute the reciprocal condition number
        double rcond = 1.0 / cond;

        // Check if rcond is NaN and replace it with a valid value
        if (std::isnan(rcond)) {
            rcond = 0.0; // Replace with your desired default value
        }

        // Check conditioning of J, if bad use pinv to solve the system.
        if (rcond < 2 * eps) {
            if (J.hasNaN()) {
                std::cout << "J contains NaN" << std::endl;
            }
            if (phi.hasNaN()) {
                std::cout << "phi contains NaN" << std::endl;
            }

            Eigen::MatrixXd pinv_J;

            try {
                // Use pinv to solve the system
                pinv_J = J.completeOrthogonalDecomposition().pseudoInverse();
            } catch (const std::exception& e) {
                std::cout << "Exception caught: " << e.what() << std::endl;
            }

            Eigen::VectorXd dx = pinv_J * (-phi);
            assert(dx.size() == N && dx.cols() == 1 && "dx is not a column vector");

            assert(dx.allFinite() && "dx is not real");
        } else {
            Eigen::SparseMatrix<double> sparse_J = J.sparseView();
            Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
            solver.analyzePattern(sparse_J);
            solver.factorize(sparse_J);

            if (solver.info() != Eigen::Success) {
                std::cout << "Failed to factorize the matrix J" << std::endl;
                // Handle the error accordingly
            }

            Eigen::VectorXd dx = solver.solve(-phi);
            assert(dx.size() == N && dx.cols() == 1 && "dx is not a column vector");
            assert(dx.allFinite() && "dx is not real");

            Eigen::MatrixXd dense_J = sparse_J;
            J = dense_J;
        }
        //! END PERTURBATION

        // tests whether the search direction is below machine precision
        auto nabla_phi = phi.transpose()*J;
        assert(nabla_phi.size() == N);

        if (dx.lpNorm<Eigen::Infinity>() < eps)
        {
            flag = 5;
            break;
        }

        // test whether we are stuck in a local minima
        if (nabla_phi.norm() < tol_abs)
        {
            flag = 6;
            break;
        }

        // test whether our direction is a sufficient descent direction
        if ((nabla_phi.dot(dx) > -rho*dx.dot(dx)))
        {
            flag = 100;
            break;
        }

        searchstrategymonitor.reset();

        auto f_0 = err;
        auto grad_f = nabla_phi.dot(dx);
        searchstrategymonitor.init( f_0, grad_f );
        searchstrategymonitor.reset();
        VectorXd x_k = x;
        assert(x_k.size() == N);
        double f_k = f_0;
        while(!searchstrategymonitor.finished(f_k))
        {
            searchstrategymonitor.step();
            auto tau = searchstrategymonitor.step_length();
            x_k = x.cwiseMax(VectorXd::Zero(x.size())) + dx * tau;
            assert(x_k.size() == N);
            auto y_k = A*x_k + b;
            assert(y_k.size() == N);
            auto phi_k = fischer(y_k,x_k);
            assert(phi_k.size() == N);
            f_k = compute_objective(phi_k);
        }
        x = x_k;
        assert(x.size() == 1);
        iterate++;
    }
    if(iterate >= max_iter)
    {
        iterate -= 1;
        flag = 8;
    }
}

MatrixXd getDirectionVectors(const Vector3d& N) {
    MatrixXd D(3, 4);

    // Normalize the normal vector
    Vector3d normalizedN = N.normalized();

    // random
    Vector3d random_vector = Vector3d::Random();
    random_vector.normalize();

    // first random vector
    Vector3d V1 = normalizedN.cross(random_vector).normalized();
    assert(abs(V1.dot(normalizedN)) < 1e-9);

    // Compute the second direction vector using cross product
    Vector3d V2 = normalizedN.cross(V1).normalized();
    assert(abs(V2.dot(normalizedN)) < 1e-9);

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

bool LCP::setup_lcp(VectorXd& tau, vector<Contact*>& collisions)
{

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
    MatrixXd N(numDOF,numCollisions); N.setZero();
    MatrixXd B(numDOF,basis_dim*numCollisions); B.setZero() ; // four basis vectors for friction
    MatrixXd E(numCollisions*basis_dim,numCollisions); E.setZero();
    MatrixXd Mu(numCollisions,numCollisions); Mu.setZero();
    Mu.setIdentity(); Mu *= mu_;

    // friction cone basis matrix
    for (int colIdx = 0; colIdx < numCollisions; colIdx ++)
    {
        Vector3d normal = collisions[colIdx]->n;
        MatrixXd D = getDirectionVectors(normal);

        //
        int objectIdxA = collisions[colIdx]->objectIdxA;
        int objectIdxB = collisions[colIdx]->objectIdxB;

        //
        MatrixXd Ja; MatrixXd Jb;

        // contact jacobians
        if (collisions[colIdx]->objectIdxA == 0)
        {
            Ja = -contactJacobian(normal, Vector3d(0, 0, 0));
            Jb = contactJacobian(normal, Vector3d(0, 0, 0));
        }
        else if (collisions[colIdx]->objectIdxA == 1)
        {
            Ja = contactJacobian(normal, Vector3d(0, 0, 0));
            Jb = -contactJacobian(normal, Vector3d(0, 0, 0));
        }

        N.block<3, 1>(3 * objectIdxA, colIdx) = Ja * normal;
        N.block<3, 1>(3 * objectIdxB, colIdx) = Jb * normal;

        // add to tangential jacobian
        B.block<3,4>(3*objectIdxA,basis_dim*colIdx) = Ja*D;
        B.block<3,4>(3*objectIdxB,basis_dim*colIdx) = Jb*D;

        // fill e matrix
        E.block<4,1>(basis_dim*colIdx,colIdx) = Vector4d(1,1,1,1);
    }

    // A Matrix
    A.resize(numCollisions + numCollisions*basis_dim + numCollisions,numCollisions + numCollisions*basis_dim + numCollisions);
    A.setZero();
    A.block(0,0,numCollisions,numCollisions) = dt*N.transpose()*M_inv*N;
    A.block(0,numCollisions,numCollisions,numCollisions*basis_dim) = dt*N.transpose()*M_inv*B;
    A.block(0,numCollisions*basis_dim + numCollisions,numCollisions,numCollisions) = MatrixXd(numCollisions,numCollisions).setZero();
    A.block(numCollisions,0,basis_dim*numCollisions,numCollisions) = dt*B.transpose()*M_inv*N;
    A.block(0,numCollisions,basis_dim*numCollisions,numCollisions*basis_dim) = dt*B.transpose()*M_inv*B;
    A.block(0,numCollisions,basis_dim*numCollisions,numCollisions) = E;
    A.block(numCollisions + basis_dim*numCollisions,0,numCollisions,numCollisions) = Mu;
    A.block(numCollisions + basis_dim*numCollisions,numCollisions,numCollisions,numCollisions*basis_dim) = -E.transpose();
    A.block(numCollisions + basis_dim*numCollisions,numCollisions*basis_dim + numCollisions,numCollisions,numCollisions) = MatrixXd(numCollisions,numCollisions).setZero();

    // b vector
    b.resize(numCollisions + numCollisions*basis_dim + numCollisions); B.setZero();
    b.head(numCollisions) = N.transpose()*M_inv*tau;
    b.segment(numCollisions,numCollisions*basis_dim ) = B.transpose()*M_inv*tau;
    b.tail(numCollisions) = VectorXd(numCollisions).setZero();

    return 1;
}

void LCP::solve_lcp_lemke(VectorXd &lambda)
{
    lambda.resize(b.size());
    lambda.setZero();
//    lemkeSolver(A,b,*lambda);

    VectorXd* lambdaPtr = &lambda;

    // Create an instance of the MobyLCPSolver class
    drake::solvers::MobyLCPSolver<double> lcp_solver;
    if(!lcp_solver.SolveLcpLemkeRegularized(A,b,lambdaPtr))
        cout << " Solver failed to find a solution" << endl;

    cout << lambda.transpose() << endl;
}

void LCP::solve_lcp_newton(Eigen::VectorXd &lambda)
{
    lambda.resize(b.size());
    lambda.setZero();
    fischer_newton_solver(A,b,lambda);
}