#include <Eigen/Dense>
#include <iostream>

using namespace Eigen;
using namespace std;

inline VectorXd gauss_seidel(const MatrixXd& A, const VectorXd& b, VectorXd& lambda, double tol = 1e-6, int max_iter = 10) {
    int n = A.rows();
    lambda = VectorXd::Zero(n);
    VectorXd residual = b - A * lambda;
    int iter = 0;

    while (residual.norm() > tol && iter < max_iter) {
        for (int j = 0; j < n; j++) {
            double sum_a_j = 0;
            for (int k = 0; k < n; k++) {
                if (k != j) {
                    sum_a_j += A(j, k) * lambda[k];
                }
            }
            lambda[j] = (b[j] - sum_a_j) / A(j, j);
        }
        residual = b - A * lambda;
        iter++;
    }

//    if (iter == max_iter) {
//        std::cout << "Gauss-Seidel reached the maximum number of iterations without converging." << std::endl;
//    } else {
//        std::cout << "Gauss-Seidel converged in " << iter << " iterations." << std::endl;
//    }

    return residual;
}



inline void jacobi(MatrixXd A, VectorXd b, VectorXd &lambda, double tol = 1e-6, int max_iter = 10) {
    int n = A.rows();
    lambda = VectorXd::Zero(n);

    for (int i = 0; i < max_iter; i++) {
        VectorXd lambda_new = VectorXd::Zero(n);
        for (int j = 0; j < n; j++) {
            double sum_a_j = 0;
            for (int k = 0; k < n; k++) {
                if (k != j) {
                    sum_a_j += A(j, k) * lambda[k];
                }
            }
            lambda_new[j] = (b[j] - sum_a_j) / A(j, j);
        }

        if ((lambda - lambda_new).norm() < tol)
        {
            break;
        }

        lambda = lambda_new;
    }
}
