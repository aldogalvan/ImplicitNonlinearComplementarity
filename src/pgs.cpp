
#include "pgs.hpp"


void ProjectedGaussSeidel::solve_projected_gauss_seidel(Eigen::VectorXd &q_n, const Eigen::VectorXd &q_n_unconstrained,
                                                        Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                        vector<Contact *> &contacts, double dt, int max_it,
                                                        double tol)
{
    // the initial states of the system
    VectorXd p = q_n_unconstrained.block<3,1>(0,0);
    Quaterniond q = Quaterniond (q_n_unconstrained(3),q_n_unconstrained(4),
                                 q_n_unconstrained(5),q_n_unconstrained(6));
    VectorXd v = u_tilde.block<3,1>(0,0);
    VectorXd w = u_tilde.block<3,1>(3,0);

    VectorXd lambda(contacts.size()); lambda.setZero();
    VectorXd delta_lambda(contacts.size()); delta_lambda.setZero();

    for (int n = 0; n < max_it; n++)
    {
        for (int cidx = 0; cidx < contacts.size(); cidx++)
        {
            
        }
    }
}

void ProjectedGaussSeidel::solve_projected_gauss_seidel_rigid(Eigen::VectorXd &q_n, const Eigen::VectorXd &q_n_unconstrained,
                                                        Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                        vector<Contact *> &contacts, double dt, int max_it,
                                                        double tol)
{
    for (int n = 0; n < max_it; n++)
    {
        for (int cidx = 0; cidx < contacts.size(); cidx++)
        {

        }
    }
}