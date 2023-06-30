#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_IMPLICIT_LCP_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_IMPLICIT_LCP_HPP

#include <Eigen/Dense>
#include <vector>
#include "chai3d.h"
#include  "objects.hpp"
#include "collision.hpp"

using namespace Eigen;
using namespace std;
using namespace chai3d;


class ImplicitLCP
{

public:

    ImplicitLCP(vector<RigidObject*> objects){
        m_objects = objects;
        m_collision_detector = new CollisionDetector(objects);
    }
    ImplicitLCP(){}

    ~ImplicitLCP(){}

public:

    // performs an explicit step
    void explicit_step(double dt);

    // steps the simulation forward by dt
    void step(double dt);

    // steps the simulation forward by dt
    void step_sparse(double dt);

    // steps the simulation forward by dt
    void step_matrix_free(double dt);

    // this function will compute
    void compute_constraints(double dt);

    // this function sets up the implicit lcp problem as per
    // "Non-Smooth Newton Methods for Deformable Multi-Body Dynamics"
    // by Macklin et al.

    static void setup_implicit_lcp_rigid(VectorXd& q_n, const VectorXd& q_n_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde, const VectorXd& M,
                                         vector<Contact*> &contacts, double dt = 0.001, int newton_it = 20, int linear_it = 20,
                                         double newton_tol = 1e-6, double linear_tol = 1e-6);

    static void setup_implicit_lcp_deformable(MatrixXd& q_n, const MatrixXd& q_unconstrained, MatrixXd &u_n, const MatrixXd &u_tilde,
                                                    const VectorXd& M, const MatrixXi& verts, const MatrixXi& tets,
                                                    vector<Contact*> &contacts, double dt, int newton_it, int linear_it,
                                                    double newton_tol, double linear_tol);

    static void setup_implicit_lcp_haptic_quasistatic(VectorXd& q_n, const Eigen::VectorXd &q_tilde, Eigen::VectorXd &u_n,
                                   vector<Contact*> &contacts, double dt = 0.001, int newton_it = 5, int linear_it = 10,
                                   double newton_tol = 1e-6, double linear_tol = 1e-6);

    static void setup_implicit_lcp_haptic_no_friction(VectorXd& q_n, const VectorXd& q_n_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                      vector<Contact*> &contacts, double dt = 0.001, int newton_it = 20, int linear_it = 20,
                                                      double newton_tol = 1e-6, double linear_tol = 1e-6);

    static void setup_implicit_lcp_haptic_friction(VectorXd& q_n, const VectorXd& q_n_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                      vector<Contact*> &contacts, double dt = 0.001, int newton_it = 20, int linear_it = 20,
                                                      double newton_tol = 1e-6, double linear_tol = 1e-6);



    //PUBLIC MEMBERS
public:

    // the objects in this simulation
    vector<RigidObject*> m_objects;

    // Matrix Assembly
    MatrixXd M;
    MatrixXd H;
    MatrixXd J;
    MatrixXd C;
    VectorXd g;
    VectorXd h;

    // iterations
    int newton_iterations = 10;
    double newton_tol = 1e-9;
    int linear_iterations = 10;
    double linear_tol = 1e-9;
    VectorXd r_normal , r_friction; //preconditioning parameter

    // the collision detector
    CollisionDetector* m_collision_detector;

};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_IMPLICIT_LCP_PPH
