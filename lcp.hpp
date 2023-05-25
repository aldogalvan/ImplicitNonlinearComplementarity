#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP

#include <Eigen/Dense>
#include "collision.hpp"
#include "objects.hpp"


class LCP
{
public:

    LCP(vector<RigidObject*> bodies)
    {
        rigid_bodies = bodies;
        M_inv = MatrixXd(bodies.size(),bodies.size()).setZero();
        M_inv = M_inv.inverse();
    }

    ~LCP(){}

    bool setup_lcp(VectorXd& tau, vector<ColInfo*>& collisions);
    void solve_lcp_newton(VectorXd& lambda);
    void solve_minimum_map(VectorXd& lambda);
    void solve_lcp_lemke(VectorXd* lambda);

    // the rigid bodies
    vector<RigidObject*> rigid_bodies;

    // Inverse Mass matrix
    MatrixXd M_inv;
    double mu_;

    //
    MatrixXd A;
    VectorXd b;
};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
