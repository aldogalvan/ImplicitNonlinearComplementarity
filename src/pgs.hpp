
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_PGS_H
#define IMPLICITNONLINEARCOMPLEMENTARITY_PGS_H

#include <Eigen/Dense>
#include <vector>
#include "chai3d.h"
#include  "objects.hpp"
#include "collision.hpp"

class ProjectedGaussSeidel
{
public:

    ProjectedGaussSeidel(){};
    ~ProjectedGaussSeidel(){};

    static void solve_projected_gauss_seidel(VectorXd& q_n, const VectorXd& q_n_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                   vector<Contact*> &contacts, double dt = 0.001, int max_it = 100, double tol = 1e-6);

    static void solve_projected_gauss_seidel_rigid(Eigen::VectorXd &q_n, const Eigen::VectorXd &q_n_unconstrained,
                                                                  Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                                  vector<Contact *> &contacts, double dt, int max_it,
                                                                  double tol);
};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_PGS_H
