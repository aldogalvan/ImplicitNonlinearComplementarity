#include "implicit_lcp.hpp"
#include "solvers.cpp"
#include "helper.hpp"
#include <iostream>
#include <cmath>


inline MatrixXd tangent_basis(const Vector3d& normal)
{
    MatrixXd D(3, 2);

    // Normalize the normal vector
    Vector3d normalizedN = normal.normalized();

    // random
    Vector3d random_vector = Vector3d::Random();
    random_vector.normalize();

    // first random vector
    Vector3d V1 = normalizedN.cross(random_vector).normalized();
    assert(abs(V1.dot(normalizedN)) < 1e-12);

    // Compute the second direction vector using cross product
    Vector3d V2 = normalizedN.cross(V1).normalized();

    assert(abs(V2.dot(normalizedN)) < 1e-12);

    // Store the direction vectors in the matrix D
    D.col(0) = V1;
    D.col(1) = V2;

    return D;
}




inline double W_fischer(const MatrixXd& D, const VectorXd& qdot, const VectorXd& lambda_f,
                   const double& r, const double& mu, const double& lambda_n)
{
    double ret = 1;
    if (lambda_n * mu > 0)
    {
        ret *= r*(sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)) - r*(mu*lambda_n - lambda_f.norm()))/
                ((D.transpose()*qdot).squaredNorm()+r*mu*lambda_n-sqrt((D.transpose()*qdot).squaredNorm()+pow(r,2)*pow((mu*lambda_n - lambda_f.norm()),2)));
    }
    return ret;
}

inline double fischer(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    return a + b - sqrt(pow(a,2) + pow(b,2)) ;
}

inline double fischer_partial_a(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    double ret;
    if (abs(a) < 1e-12 && abs(b) < 1e-12)
    {
        ret = 0;
    }
    else
    {
        ret = 1 - (a/sqrt(pow(a,2) + pow(b,2)));
    }
    return ret;
}

inline double fischer_partial_b(const double& a, const double& b)
{
    assert(!isnan(a) && !isnan(b));
    double ret;
    if (abs(a) < 1e-12 && abs(b) < 1e-12)
    {
        ret = 1;
    }
    else
    {
        ret = 1 - (b/sqrt(pow(a,2) + pow(b,2)));
    }
    return ret;
}

void ImplicitLCP::compute_constraints(double dt)
{
//    int num_contacts = m_collision_detector->m_contacts.size();
//
//    ////////////////////////////////////////////////////////////////////
//    //////// IF WE HAVE COLLISIONS/CONSTRAINTS /////////////////////////
//    ////////////////////////////////////////////////////////////////////
//
//    if (num_contacts > 0) {
//
//        int num_objects = m_objects.size();
//
//        MatrixXd K = MatrixXd(6*num_objects,6*num_objects); K.setZero();
//        VectorXd phi_normal(num_contacts); phi_normal.setZero();
//        VectorXd phi_friction(2 * num_contacts); phi_friction.setZero();
//        VectorXd lambda_normal(num_contacts); lambda_normal.setZero();
//        VectorXd lambda_friction(2 * num_contacts); lambda_friction.setZero();
//        VectorXd lambda(lambda_normal.size()); lambda.setZero();
////        VectorXd lambda(lambda_normal.size() + lambda_friction.size()); lambda.setZero();
//        MatrixXd J_normal(num_contacts, num_objects * 6); J_normal.setZero();
//        MatrixXd J_friction(2 * num_contacts, num_objects * 6); J_friction.setZero();
//        MatrixXd S(num_contacts, num_contacts); S.setZero();
//        MatrixXd W(2 * num_contacts, 2 * num_contacts); W.setZero();
//        VectorXi is_active(num_contacts); // determines if constraint is active or not
//
//        // the mass matrix
//        MatrixXd M_tilde(6 * num_objects, 6 * num_objects);
//        M_tilde.setZero();
//        MatrixXd M_tilde_inv(6 * num_objects, 6 * num_objects);
//        M_tilde_inv.setZero();
//
//        VectorXd q(7 * num_objects); // the constrained position in global coordinates
//        VectorXd u(6 * num_objects); // the constrained velocity in generalized coordinates
//
//        for (int bidx = 0; bidx < num_objects; bidx++)
//        {
//            MatrixXd G = m_objects[bidx]->kinematic_map_G();
//            M_tilde.block<6, 6>(6 * bidx, 6 * bidx) = m_objects[bidx]->mass_matrix();
//            q.block<3,1>(6*bidx,0) = m_objects[bidx]->x_unconstrained;
//            q(6*bidx+3) = m_objects[bidx]->q_unconstrained.w();
//            q(6*bidx+4) = m_objects[bidx]->q_unconstrained.x();
//            q(6*bidx+5) = m_objects[bidx]->q_unconstrained.y();
//            q(6*bidx+6) = m_objects[bidx]->q_unconstrained.z();
//            u.block<3,1>(6*bidx,0) = m_objects[bidx]->xdot;
//            u.block<3,1>(6*bidx + 3,0) = m_objects[bidx]->omega;
//        }
//
//        VectorXd u_tilde = u;
//
//        ////////////////////////////////////////////////////////////////////////
//        ////////////////// START THE NEWTON ITERATIONS ////////////////////////
//        ////////////////////////////////////////////////////////////////////////
//
//        for (int n = 0; n < newton_iterations; n++)
//        {
//            ////////////////////////////////////////////////////////////////////////
//            ////////////////// ASSEMBLE THE CONTACTS ///////////////////////////////\
//            ///////////////////////////////////////////////////////////////////////
//
//            for (int cidx = 0; cidx < num_contacts; cidx++) {
//
//                // the collision
//                auto &contact = m_collision_detector->m_contacts[cidx];
//
//                // object indices
//                auto &objectA = m_objects[contact->bodyIdxA];
//                auto &objectB = m_objects[contact->bodyIdxB];
//
//                // set the contact index
//                int body_idxA = contact->bodyIdxA;
//                int body_idxB = contact->bodyIdxB;
//
//                // velocities and positins
//                VectorXd u_A = u.block<6,1>(6*body_idxA,0);
//                VectorXd q_A = q.block<7,1>(7*body_idxA,0);
//                VectorXd u_B = u.block<6,1>(6*body_idxB,0);
//                VectorXd q_B = q.block<7,1>(7*body_idxB,0);
//
//                // set the preconditioner values
//                double r_friction, r_normal = 1;
//
//                // collision values
//                Vector3d normal = contact->normal;
//                cout << "DEBUG : normal = " << normal.transpose() << endl;
//                Vector3d contact_wrt_bodyA = contact->contact_wrt_bodyA;
//                Vector3d contact_wrt_bodyB = contact->contact_wrt_bodyB;
//
//
//                // body A configuration
//                Vector3d pos_A = q_A.head(3);
//                Quaterniond rot_A(q_A(3),q_A(4),q_A(5),q_A(6));
//                rot_A.normalize();
//
//                // body B configuration
//                Vector3d pos_B = q_B.head(3);
//                Quaterniond rot_B(q_B(3),q_B(4),q_B(5),q_B(6));
//                rot_B.normalize();
//
//                // get the contact pts
//                Vector3d contact_ptA = rot_A * contact_wrt_bodyA + pos_A;
//                cout << "DEBUG: contact_wrt_bodyA = " << contact_wrt_bodyA.transpose() << endl;
//                cout << "DEBUG: contact_ptA = " << contact_ptA.transpose() << endl;
//                Vector3d contact_ptB = rot_B * contact_wrt_bodyB + pos_B;
//                cout << "DEBUG: contact_wrt_bodyB = " << contact_wrt_bodyB.transpose() << endl;
//                cout << "DEBUG: contact_ptB = " << contact_ptB.transpose() << endl;
//
//
//                // compute penetration depth (should already be done but w/e)
//                double penetration_depth = (contact_ptA - contact_ptB).dot(normal);
//                cout << "DEBUG : penetration_depth = " << penetration_depth << endl;
//
//                // if the penetration depth is > 0 constraint is inactive
//                if (penetration_depth > 0) {
//                    is_active[cidx] = 0;
//                } else {
//                    is_active[cidx] = 1;
//                }
//
//                // the kinematic map matrices
//                MatrixXd G_A = m_objects[body_idxA]->kinematic_map_G();
//                cout << "DEBUG: G_A = " << G_A << endl;
//                MatrixXd G_B = m_objects[body_idxB]->kinematic_map_G();
//                cout << "DEBUG: G_B = " << G_B << endl;
//
//                // the constraint value
//                double C_n = penetration_depth;
//
//                // gradients of the constraint for normal
//                VectorXd grad_C_n_A(6); grad_C_n_A.head(3) = normal; grad_C_n_A.tail(3) = vec2skew(-normal).transpose()*contact_wrt_bodyA;
//                cout << "DEBUG: grad_C_n_A  = " << grad_C_n_A.transpose() << endl;
//                VectorXd grad_C_n_B(6); grad_C_n_B.head(3) = -normal; grad_C_n_B.tail(3) = vec2skew(normal).transpose()*contact_wrt_bodyB;
//                cout << "DEBUG: grad_C_n_B  = " << grad_C_n_B.transpose() << endl;
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Normal Constraint Matrices for Object A ///////
//                //////////////////////////////////////////////////////////////////////
//
//                VectorXd phi_normal_wrt_q = fischer_partial_a(C_n, r_normal * lambda_normal(cidx)) * grad_C_n_A;
//                cout << "DEBUG: " << phi_normal_wrt_q << endl;
//                assert(!phi_normal_wrt_q.hasNaN());
//                double phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal * lambda_normal(cidx)) * r_normal;
//                assert(!isnan(phi_normal_wrt_lambda));
//                S(cidx, cidx) = phi_normal_wrt_lambda;
//                J_normal.block<1, 6>(cidx, 6*body_idxA) = phi_normal_wrt_q.transpose();
//
//                // fill the normal constraint vector and normal jacobian
//                phi_normal(cidx) = fischer(C_n, r_normal * lambda_normal(cidx));
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Normal Constraint Matrices for Object B ///////
//                //////////////////////////////////////////////////////////////////////
//
//                phi_normal_wrt_q = fischer_partial_a(C_n, r_normal * lambda_normal(cidx)) * grad_C_n_B;
//                assert(!phi_normal_wrt_q.hasNaN());
//                phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal * lambda_normal(cidx)) * r_normal;
//                assert(!isnan(phi_normal_wrt_lambda));
//                S(cidx, cidx) = phi_normal_wrt_lambda;
//                J_normal.block<1, 6>(cidx, 6*body_idxB) = phi_normal_wrt_q.transpose();
//
//                // fill the normal constraint vector and normal jacobian
//                phi_normal(cidx) = fischer(C_n, r_normal * lambda_normal(cidx));
//
//
//                // Friction basis
//                MatrixXd D_n_A = tangent_basis(normal); MatrixXd D_n_B = -D_n_A;
//                D_n_B.conservativeResize(6,2);
//                D_n_B.block<3,1>(3,0)= vec2skew(-D_n_B.block<3,1>(0,0)).transpose()*contact_wrt_bodyB;
//                D_n_B.block<3,1>(3,1)= vec2skew(-D_n_B.block<3,1>(0,1)).transpose()*contact_wrt_bodyB;
//                D_n_A.conservativeResize(6,2);
//                D_n_A.block<3,1>(3,0)= vec2skew(-D_n_A.block<3,1>(0,0)).transpose()*(contact_wrt_bodyA);
//                D_n_A.block<3,1>(3,1)= vec2skew(-D_n_A.block<3,1>(0,1)).transpose()*(contact_wrt_bodyA);
//
//                ////////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Friction Constraint Matrices for Object A ///////
//                ////////////////////////////////////////////////////////////////////////
//
//
//                auto phi_friction_wrt_qdot = D_n_A;
//                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
//                double W_f = W_fischer(D_n_A, u_A, lambda_f, r_friction, 0.1, lambda_normal(cidx));
//                MatrixXd phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();
//
//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_A.transpose() * u_A + W_f * lambda_f;
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Friction Constraint Matrices for Object B ///////
//                //////////////////////////////////////////////////////////////////////
//
//
//                phi_friction_wrt_qdot = D_n_B;
//                lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
//                W_f = W_fischer(D_n_B, u_B, lambda_f, r_friction, 0.1, lambda_normal(cidx));
//                phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();
//
//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_B.transpose() * u_B + W_f * lambda_f;
//
//            }
//
//            assert(!J_friction.hasNaN());
//            assert(!J_normal.hasNaN());
//            assert(!S.hasNaN());
//            assert(!W.hasNaN());
//            H = M_tilde - pow(dt, 2) * K;
//            assert(!H.hasNaN());
//            assert(J_normal.cols() == J_friction.cols());
//            J.resize(J_normal.rows(), J_normal.cols());
////        J.resize(J_normal.rows() + J_friction.rows(), J_normal.cols());
//            J.setZero();
//            J << J_normal;
////        J << J_normal, J_friction;
//            assert(!J.hasNaN());
////            cout << "DEBUG : J = " << J << endl;
//
//            C.resize(S.rows(), S.cols());
////        C.resize(E.rows() + S.rows() + W.rows(), E.cols() + S.cols() + W.cols());
//            C.setZero();
//            C.block(0, 0, S.rows(), S.cols()) = S / pow(dt, 2);
////        C.block(E.rows() + S.rows(), E.cols() + S.cols(), W.rows(), W.cols()) = W / dt;
//            assert(!C.hasNaN());
////            cout << "DEBUG : C = " << C << endl;
//
//            VectorXd delta_lambda(lambda_normal.size());
////        VectorXd delta_lambda(lambda_normal.size() + lambda_friction.size());
//            delta_lambda.setZero();
//            // not sure if correct
//            // delta_lambda << lambda_normal, lambda_friction;
//            assert(!lambda.hasNaN());
//            assert(J.rows() == lambda.size());
//            assert(M_tilde.cols() == u.size());
//            g = M_tilde * (u) - dt * J.transpose() * lambda;
////        cout << "DEBUG: u " << u_n.transpose() << " , u_tilde = " << u_tilde.transpose() << endl;
//            assert(!g.hasNaN());
////        cout << "DEBUG : g = " << g.transpose() << endl;
//            MatrixXd I_n(phi_normal.size(), phi_normal.size());
//            I_n.setIdentity();
//            I_n /= dt;
//            assert(!I_n.hasNaN());
//            MatrixXd I_f(phi_friction.size(), phi_friction.size());
//            I_f.setIdentity();
//            assert(!I_f.hasNaN());
//            MatrixXd I(I_n.rows(), I_n.cols());
////        MatrixXd I(I_n.rows() + I_f.rows(), I_n.cols() + I_f.cols());
//            I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
////        cout << "DEBUG : I = " << I << endl;
////        I.block(I_n.rows(), I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
//            VectorXd phi(phi_normal.size());
////        VectorXd phi(phi_normal.size() + phi_friction.size());
//            phi << phi_normal;
////        phi << phi_normal, phi_friction;
//            assert(!phi.hasNaN());
////        cout << "DEBUG: phi = " << phi.transpose() << endl;
//            h = I * phi;
////        cout << "DEBUG: h = " << h.transpose() << endl;
//
//            assert(!h.hasNaN());
//            assert(!J.hasNaN());
//            assert(J.rows() == C.rows());
//            assert(J.rows() == C.cols());
//            auto H_inv = H.inverse();
//            // zero out any nans caused by inversion
//            assert(!H_inv.hasNaN());
//            MatrixXd A = J * H_inv * J.transpose() + C;
//
//            if (is_ill_conditioned(A,1000))
//            {
//                cout << "DEBUG: Matrix A is ill conditioned!" << endl;
//            }
//            assert(!A.hasNaN());
//            assert(J.cols() == H.rows() && H.cols() == g.size());
//            assert(J.rows() == h.size());
//            VectorXd b = (J * H_inv * g - h) / dt;
//            assert(!b.hasNaN());
//
//
//            // zero out any inactive contacts
//            for (int contact = 0; contact < num_contacts; contact++)
//            {
//                if (is_active[contact] == 0 )
//                {
//                    A.row(num_contacts + 2*contact) = RowVectorXd::Zero(A.cols());
//                    A.row(num_contacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
//                    b(num_contacts + 2*contact) = 0;
//                    b(num_contacts + 2*contact + 1) = 0;
//                }
//            }
//
//            // solve this system with Gauss-Seidel
//            auto res = gauss_seidel(A, b, delta_lambda,linear_tol);
////        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;
////        cout << "DEBUG: residual = " << res.transpose();
//            assert(!delta_lambda.hasNaN());
//            auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
//            assert(!delta_u.hasNaN());
//            // perform a line search
//            int t = 1;
//
//            // update solution
//            lambda += t*delta_lambda;
//            u += t*delta_u;
//            for (int bidx = 0 ; bidx < m_objects.size(); bidx++)
//            {
//                q.block<7,1>(bidx*7,0) += dt*m_objects[bidx]->kinematic_map_G()*u.block<7,1>(bidx*7,0);
//                m_objects[bidx]->x = q.block<3,1>(bidx*7,0); m_objects[bidx]->q = q.block<4,1>(bidx*7+3,0);
//                m_objects[bidx]->xdot = u.block<3,1>(bidx*6,0); m_objects[bidx]->omega = q.block<3,1>(bidx*6+3,0);
//            }
//
//            // set the constraint forces
//            lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
////            lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);
//
//            // if the residual is low enough then break out
//            if (res.norm() < newton_tol)
//            {
////                cout << "Newtons Method Converged in " << n << " iterations!" << endl;
//                return;
//            }
//
//            r_normal = pow(dt,2)*(J*M_tilde_inv*J.transpose()).diagonal();
//            r_friction = dt*(J*M_tilde_inv*J.transpose()).diagonal();
//        }
//    }
//    else
//    {
//        for (auto b: m_objects)
//        {
//            b->x = b->x_unconstrained;
//            b->q = b->q_unconstrained;
//        }
//        return;
//    }
}

void ImplicitLCP::step(double dt)
{

//    for(auto b : m_objects) {
//        b->f = b->mass * Eigen::Vector3d(0., 0, -.0098);
//        b->tau.setZero();
//        b->fc.setZero();
//        b->tauc.setZero();
//    }
//
//    for (auto c : m_collision_detector->m_contacts)
//    {
//        delete [] c;
//    }
//    m_collision_detector->m_contacts.clear();
//
//    // update the inertia matrices
//    for(auto b : m_objects)
//    {
//        b->update_inertia_matrix();
//    }
//
//    // compute collisions
//    m_collision_detector->computeCollisions();
//
//    // compute the constraint impulses
//    compute_constraints(dt);
//
//    // time integrate the objects in the system
//    for (auto b : m_objects)
//    {
//        if (!b->is_static)
//        {
//            b->xdot += dt * (1 / b->mass) * b->f;
//            b->omega += dt * b->Iinv * (b->tau - b->omega.cross(b->I * b->omega));
//            b->x_unconstrained += dt * b->xdot;
//            VectorXd delta_q = dt * (0.5 * b->kinematic_map_H() * b->omega);
//            b->q_unconstrained.w() += Quaterniond(delta_q(0), delta_q(1), delta_q(2), delta_q(3)).w();
//            b->q_unconstrained.vec() += Quaterniond(delta_q(0), delta_q(1), delta_q(2), delta_q(3)).vec();
//            b->q_unconstrained.norm();
//        }
//    }
}



void ImplicitLCP::setup_implicit_lcp_haptic_no_friction(VectorXd& q_n, const VectorXd& q_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                        vector<Contact*> &contacts, double dt, int newton_it, int linear_it,
                                                        double newton_tol, double linear_tol)
{
    assert(!u_tilde.hasNaN()); assert(!q_n.hasNaN());

    // linear translation
    int ndof = 3;

    MatrixXd M = Matrix3d().setIdentity()*100;
    MatrixXd M_inv = M.inverse();
    int numContacts = contacts.size();
    VectorXd phi_normal(numContacts); phi_normal.setZero();
    VectorXd phi_friction(2*numContacts); phi_friction.setZero();
    VectorXd lambda_normal(numContacts); lambda_normal.setZero();
    VectorXd lambda_friction(2*numContacts); lambda_friction.setZero();
    VectorXd lambda(lambda_normal.size()); lambda.setZero();
    MatrixXd J_normal(numContacts,ndof); J_normal.setZero();
    MatrixXd J_friction(2*numContacts,ndof); J_friction.setZero();
    MatrixXd E;
    MatrixXd S(numContacts,numContacts); S.setZero();
    MatrixXd W(2*numContacts,2*numContacts); W.setZero();
    VectorXi is_active(numContacts);
    VectorXd r_normal = VectorXd::Ones(numContacts);
    VectorXd r_friction = VectorXd::Ones(numContacts);

    for (int n = 0 ; n < newton_it; n++) {

        cout << "NEWTON IT: " << n << endl;

        // assemble the matrix
        for (int colIdx = 0; colIdx < numContacts; colIdx++) {

            // set the preconditioner values
            double r_friction_ = dt; double r_normal_ = dt*dt;

            // collision values
            Vector3d normal = contacts[colIdx]->n;
            Vector3d contact_wrt_bodyA = contacts[colIdx]->contact_wrt_bodyA;
            Vector3d contact_wrt_bodyB = contacts[colIdx]->contact_wrt_bodyB;

            // set the contact index
            int contact_idx_bodyA = contacts[colIdx]->bodyIdxA;
            int contact_idx_bodyB = contacts[colIdx]->bodyIdxB;

            // get the contact pts
            //! need to do the full rotation transform
            Vector3d contact_ptA;
            if (contact_idx_bodyA != -1)
            {
                contact_ptA = q_n.block<3, 1>(contact_idx_bodyA, 0) + contact_wrt_bodyA;
            }
            else
            {
                contact_ptA =  contact_wrt_bodyA;
            }
            Vector3d contact_ptB;
            if (contact_idx_bodyB != -1 )
            {
                contact_ptB = q_n.block<3,1>(contact_idx_bodyB,0) + contact_wrt_bodyB;
            }
            else{
                contact_ptB =  contact_wrt_bodyB;
            }

            Vector3d rel_vel = u_n;
//            double penetration_depth = ().dot(normal) - 1e-6;
//            double penetration_depth = (rel_vel).dot(normal) - 1e-6;
//
//            if (penetration_depth > 0 )
//            {
//                is_active[colIdx] = 0;
//            }
//            else
//            {
//                is_active[colIdx] = 1;
//            }

            cout << "DEBUG: " << rel_vel.transpose() << endl;

            VectorXd grad_C_n_A = normal.normalized();
            VectorXd grad_C_n_B = -normal.normalized();

            //
            if (contact_idx_bodyA == 0)
            {
                double C_n = rel_vel.dot(grad_C_n_A);

                cout << "DEBUG: C_n = " << C_n << endl;

                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_A;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = phi_normal_wrt_q.transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));
            }
            if (contact_idx_bodyB == 0)
            {
                double C_n = rel_vel.dot(grad_C_n_B);

                cout << "DEBUG: C_n = " << C_n << endl;

                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_B;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = phi_normal_wrt_q.transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));

            }
        }

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = M;
        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());
        MatrixXd J(J_normal.rows(), J_normal.cols());
        J.setZero();
        J << J_normal;
        assert(!J.hasNaN());

        MatrixXd C(E.rows() + S.rows(), E.cols() + S.cols());
        C.setZero();
        C.block(0, 0, E.rows(), E.cols()) = E / pow(dt, 2);
        C.block(E.rows(), E.cols(), S.rows(), S.cols()) = S / pow(dt, 2);
        assert(!C.hasNaN());

        VectorXd delta_lambda(lambda_normal.size());
        delta_lambda.setZero();
        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(M.cols() == u_n.size());
        VectorXd g =  M*(u_n - u_tilde) - dt * J.transpose() * lambda;
        assert(!g.hasNaN());
        MatrixXd I_n(phi_normal.size(), phi_normal.size());
        I_n.setIdentity();
        I_n /= dt;
        assert(!I_n.hasNaN());
        MatrixXd I_f(phi_friction.size(), phi_friction.size());
        I_f.setIdentity();
        assert(!I_f.hasNaN());
        MatrixXd I(I_n.rows(), I_n.cols());
        I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
        VectorXd phi(phi_normal.size());
        phi << phi_normal;
        assert(!phi.hasNaN());
        VectorXd h = I * phi;

        assert(!h.hasNaN());
        assert(!J.hasNaN());
        assert(J.rows() == C.rows());
        assert(J.rows() == C.cols());
        auto H_inv = H.inverse();
        assert(!H_inv.hasNaN());
        MatrixXd A = J * H_inv * J.transpose() + C;
        A += MatrixXd::Identity(A.rows(),A.cols())*1e-6;
        if (is_ill_conditioned(A,1000))
        {
            cout << "DEBUG: Matrix A is ill conditioned!" << endl;
        }
        assert(!A.hasNaN());
        assert(J.cols() == H.rows() && H.cols() == g.size());
        assert(J.rows() == h.size());
        VectorXd b = (J * H_inv * g - h) / dt;
//        cout << "DEBUG: b = " << b.transpose() << endl;
//        cout << "DEBUG: A = " << A << endl;
        assert(!b.hasNaN());
        for (int contact = 0; contact < numContacts; contact++)
        {
            if (is_active[contact] == 0 )
            {
                // zero the normal
                A.row(contact) = RowVectorXd::Zero(A.cols());
                b(contact) = 0;

                A.row(numContacts + 2*contact) = RowVectorXd::Zero(A.cols());
                A.row(numContacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
                b(numContacts + 2*contact) = 0;
                b(numContacts + 2*contact + 1) = 0;
            }
        }

//        cout << "DEBUG: is_active = " << is_active.transpose() << endl;
        // solve this system with Gauss-Seidel
        auto res = gauss_seidel(A, b, delta_lambda,linear_tol);
//        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;
//        cout << "DEBUG: convergence = " << delta_lambda.norm() << endl;
//        cout << "DEBUG: residual = " << res.transpose();
        assert(!delta_lambda.hasNaN());
        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
//        cout << "DEBUG: delta_u = " <<  delta_u.transpose() <<  endl;
        assert(!delta_u.hasNaN());
        // perform a line search
        double t = 0.1;
        // update solution
        lambda += t*delta_lambda;
//        cout << delta_lambda.transpose() << endl;
        u_n += t*delta_u;
        q_n += dt*u_n;

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);

        // if the residual is low enough then break out
        r_normal = pow(dt,2)*(J*M_inv*J.transpose()).diagonal();
        r_friction = dt*(J*M_inv*J.transpose()).diagonal();
    }

//    cout << "Newtons Method Failed to Converge" << endl;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void ImplicitLCP::setup_implicit_lcp_haptic_friction(VectorXd& q_n, const VectorXd& q_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                                        vector<Contact*> &contacts, double dt, int newton_it, int linear_it,
                                                        double newton_tol, double linear_tol)
{
    assert(!u_tilde.hasNaN()); assert(!q_n.hasNaN());

    // linear translation
    int ndof = 3;

    MatrixXd M = Matrix3d().setIdentity()*100;
    MatrixXd M_inv = M.inverse();
    int numContacts = contacts.size();
    VectorXd phi_normal(numContacts); phi_normal.setZero();
    VectorXd phi_friction(2*numContacts); phi_friction.setZero();
    VectorXd lambda_normal(numContacts); lambda_normal.setZero();
    VectorXd lambda_friction(2*numContacts); lambda_friction.setZero();
    VectorXd lambda(lambda_normal.size() + lambda_friction.size()); lambda.setZero();
    MatrixXd J_normal(numContacts,ndof); J_normal.setZero();
    MatrixXd J_friction(2*numContacts,ndof); J_friction.setZero();
    MatrixXd E;
    MatrixXd S(numContacts,numContacts); S.setZero();
    MatrixXd W(2*numContacts,2*numContacts); W.setZero();
    VectorXi is_active(numContacts);
    VectorXd r_normal = VectorXd::Ones(numContacts);
    VectorXd r_friction = VectorXd::Ones(numContacts);

//    cout << "DEBUG: u_tilde = " << u_tilde.transpose() << endl;

    for (int n = 0 ; n < newton_it; n++) {

        cout << "NEWTON IT: " << n << endl;

        // assemble the matrix
        for (int colIdx = 0; colIdx < numContacts; colIdx++) {

            // set the preconditioner values
            double r_friction_ = dt; double r_normal_ = dt*dt;

            // collision values
            Vector3d normal = contacts[colIdx]->n;
            Vector3d contact_wrt_bodyA = contacts[colIdx]->contact_wrt_bodyA;
            Vector3d contact_wrt_bodyB = contacts[colIdx]->contact_wrt_bodyB;

            // set the contact index
            int contact_idx_bodyA = contacts[colIdx]->bodyIdxA;
            int contact_idx_bodyB = contacts[colIdx]->bodyIdxB;

            // get the contact pts
            //! need to do the full rotation transform
            Vector3d contact_ptA;
            if (contact_idx_bodyA != -1)
            {
                contact_ptA = q_n.block<3, 1>(contact_idx_bodyA, 0) + contact_wrt_bodyA;
            }
            else
            {
                contact_ptA =  contact_wrt_bodyA;
            }
            Vector3d contact_ptB;
            if (contact_idx_bodyB != -1 )
            {
                contact_ptB = q_n.block<3,1>(contact_idx_bodyB,0) + contact_wrt_bodyB;
            }
            else{
                contact_ptB =  contact_wrt_bodyB;
            }

            Vector3d rel_vel = u_n;

//            double penetration_depth = (contact_ptA - contact_ptB).dot(normal) - 1e-6;
//            cout << "DEBUG: penetration_depth = " << penetration_depth << endl;
//            if (penetration_depth > 0 )
//            {
//                is_active[colIdx] = 0;
//            }
//            else
//            {
//                is_active[colIdx] = 1;
//            }

//            double C_n = penetration_depth;

            VectorXd grad_C_n_A = normal.normalized();
            VectorXd grad_C_n_B = -normal.normalized();

            // now the same for the friction jacobian
            auto D = tangent_basis(normal);

            //
            if (contact_idx_bodyA == 0)
            {
                double C_n = rel_vel.dot(grad_C_n_A);

                if (C_n > 0 )
                {
                    is_active[colIdx] = 0;
                }
                else
                {
                    is_active[colIdx] = 1;
                }

                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_A;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = phi_normal_wrt_q.transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));

                // now the same for the friction jacobian
                auto D_n_A = D;
                auto phi_friction_wrt_qdot = D_n_A;
                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * colIdx, 0);
                double W_f = W_fischer(D_n_A, u_n, lambda_f, r_friction_, 0.1, lambda_normal(colIdx));
                auto phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
                W.block<2, 2>(2 * colIdx, 2 * colIdx) = phi_friction_wrt_lambda_f;
                J_friction.block<2, 3>(2 * colIdx, 0) = phi_friction_wrt_qdot.transpose() * G;

                // fill the friction vector
                phi_friction.block<2, 1>(2 * colIdx, 0) = D.transpose() * u_n + W_f * lambda_f;
            }
            if (contact_idx_bodyB == 0)
            {
                double C_n = rel_vel.dot(grad_C_n_B);

                if (C_n > 0 )
                {
                    is_active[colIdx] = 0;
                }
                else
                {
                    is_active[colIdx] = 1;
                }

                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_B;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = phi_normal_wrt_q.transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));

                // now the same for the friction jacobian
                auto D_n_B = -D;
                auto phi_friction_wrt_qdot = D_n_B;
                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * colIdx, 0);
                double W_f = W_fischer(D_n_B, u_n, lambda_f, r_friction_, 0.1, lambda_normal(colIdx));
                auto phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
                W.block<2, 2>(2 * colIdx, 2 * colIdx) = phi_friction_wrt_lambda_f;
                J_friction.block<2, 3>(2 * colIdx, 0) = phi_friction_wrt_qdot.transpose() * G;

                // fill the friction vector
                phi_friction.block<2, 1>(2 * colIdx, 0) = D.transpose() * u_n + W_f * lambda_f;
            }
        }

        cout << "DEBUG: is_active = " << is_active.transpose() << endl;

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = M;

        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());
        MatrixXd J(J_normal.rows() + J_friction.rows(), J_normal.cols());
        J.setZero();
        J << J_normal, J_friction;
        assert(!J.hasNaN());

        MatrixXd C(E.rows() + S.rows() + W.rows(), E.cols() + S.cols() + W.cols());
        C.setZero();
        C.block(0, 0, E.rows(), E.cols()) = E / pow(dt, 2);
        C.block(E.rows(), E.cols(), S.rows(), S.cols()) = S / pow(dt, 2);
        C.block(E.rows() + S.rows(), E.cols() + S.cols(), W.rows(), W.cols()) = W / dt;
        assert(!C.hasNaN());

        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(M.cols() == u_n.size());
        VectorXd g = M * (u_n - u_tilde) - dt * J.transpose() * lambda;
        assert(!g.hasNaN());
        MatrixXd I_n(phi_normal.size(), phi_normal.size());
        I_n.setIdentity();
        I_n /= dt;
        assert(!I_n.hasNaN());
        MatrixXd I_f(phi_friction.size(), phi_friction.size());
        I_f.setIdentity();
        assert(!I_f.hasNaN());
        MatrixXd I(I_n.rows() + I_f.rows(), I_n.cols() + I_f.cols());
        I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
//        cout << "DEBUG : I = " << I << endl;
        I.block(I_n.rows(), I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
        VectorXd phi(phi_normal.size() + phi_friction.size());
        phi << phi_normal, phi_friction;
        assert(!phi.hasNaN());
//        cout << "DEBUG: phi = " << phi.transpose() << endl;
        VectorXd h = I * phi;
//        cout << "DEBUG: h = " << h.transpose() << endl;

        assert(!h.hasNaN());
        assert(!J.hasNaN());
        assert(J.rows() == C.rows());
        assert(J.rows() == C.cols());
        auto H_inv = H.inverse();
        assert(!H_inv.hasNaN());
        MatrixXd A = J * H_inv * J.transpose() + C;
        A += MatrixXd::Identity(A.rows(),A.cols())*1e-6;

        assert(!A.hasNaN());
        assert(J.cols() == H.rows() && H.cols() == g.size());
        assert(J.rows() == h.size());
        VectorXd b = (J * H_inv * g - h) / dt;
//        cout << "DEBUG: b = " << b.transpose() << endl;
//        cout << "DEBUG: A = " << A << endl;
        assert(!b.hasNaN());
        for (int contact = 0; contact < numContacts; contact++)
        {
            if (is_active[contact] == 0 )
            {
                // zero the normal
                A.row(contact) = RowVectorXd::Zero(A.cols());
                b(contact) = 0;

                // zero the friction
                A.row(numContacts + 2*contact) = RowVectorXd::Zero(A.cols());
                A.row(numContacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
                b(numContacts + 2*contact) = 0;
                b(numContacts + 2*contact + 1) = 0;
            }
        }

//        cout << "DEBUG: is_active = " << is_active.transpose() << endl;
        // solve this system with Gauss-Seidel
        VectorXd delta_lambda = A.fullPivLu().solve(b);
//        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;
//        cout << "DEBUG: convergence = " << delta_lambda.norm() << endl;
//        cout << "DEBUG: residual = " << res.transpose();
//        assert(!delta_lambda.hasNaN());

        if (delta_lambda.hasNaN())
        {
            cout << "DEBUG: A = "<< A << endl;
            cout << "DEBUG: b = " << b.transpose() << endl;
        }
        auto delta_u_normal = H_inv*(J.block(0,0,lambda_normal.size(),3).transpose()*delta_lambda.block(0,0,lambda_normal.size(),1)*dt - g);
        auto delta_u_friction = H_inv*(J.block(lambda_normal.size(),0,lambda_friction.size(),3).transpose()*delta_lambda.block(lambda_normal.size(),0,lambda_friction.size(),1)*dt - g);
        auto delta_u = delta_u_normal + delta_u_friction;
        cout << "DEBUG: delta_u_normal = " <<  delta_u_normal.transpose() <<  endl;
        cout << "DEBUG: delta_u_friction = " <<  delta_u_friction.transpose() <<  endl;
        cout << "DEBUG: dot = " << delta_u_normal.dot(delta_u_friction) << endl;

        assert(!delta_u.hasNaN());
        // perform a line search
        double t = 0.1;

        // update solution
        lambda += t*delta_lambda;
        u_n += t*delta_u;
        q_n += dt*u_n;

        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;
        cout << "DEBUG: q_n = " << q_n.transpose() << endl;
        cout << "DEBUG: u_n = " << u_n.transpose() << endl;

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
        lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);

        // if the residual is low enough then break out
        r_normal = pow(dt,2)*(J*M_inv*J.transpose()).diagonal();
        r_friction = dt*(J*M_inv*J.transpose()).diagonal();
    }

//    cout << "Newtons Method Failed to Converge" << endl;
}

void ImplicitLCP::setup_implicit_lcp_rigid(VectorXd& q_n, const VectorXd& q_unconstrained, Eigen::VectorXd &u_n, const Eigen::VectorXd &u_tilde,
                                           const VectorXd& M, vector<Contact*> &contacts, double dt, int newton_it, int linear_it,
                                           double newton_tol, double linear_tol)
{
    assert(!u_tilde.hasNaN()); assert(!q_n.hasNaN());
    int numContacts = contacts.size();
    VectorXd phi_normal(numContacts); phi_normal.setZero();
    VectorXd phi_friction(2*numContacts); phi_friction.setZero();
    VectorXd lambda_normal(numContacts); lambda_normal.setZero();
    VectorXd lambda_friction(2*numContacts); lambda_friction.setZero();
    VectorXd lambda(lambda_normal.size() + lambda_friction.size()); lambda.setZero();
    MatrixXd J_normal(numContacts,6); J_normal.setZero();
    MatrixXd J_friction(2*numContacts, 6); J_friction.setZero();
    MatrixXd S(numContacts,numContacts); S.setZero();
    MatrixXd W(2*numContacts,2*numContacts); W.setZero();
    VectorXi is_active(numContacts);
    VectorXd r_normal = VectorXd::Ones(numContacts);
    VectorXd r_friction = VectorXd::Ones(numContacts);

    // create the mass
    MatrixXd M_ = M.asDiagonal();
    MatrixXd M_inv_ = M_.inverse();

    for (int n = 0 ; n < newton_it; n++) {

        MatrixXd M_tilde = M_;
        Quaterniond qd(q_n(3),q_n(4),q_n(5),q_n(6));
        qd.normalize();
        MatrixXd G = RigidObject::kinematic_map_G(qd);
        cout << "DEBUG: G = " << G << endl;
        M_tilde.block<3,3>(3,3) = qd*M_tilde.block<3,3>(3,3)*qd.inverse();
        MatrixXd M_tilde_inv = M_inv_;
        M_tilde_inv.block<3,3>(3,3) = qd*M_tilde_inv.block<3,3>(3,3)*qd.inverse();

        // assemble the matrix
        for (int cidx = 0; cidx < numContacts; cidx++) {

            // set the preconditioner values
            double r_friction_ = dt; double r_normal_ = dt*dt;

            // the collision
            auto &contact = contacts[cidx];

            // set the contact index
            int body_idxA = contact->bodyIdxA;
            int body_idxB = contact->bodyIdxB;

//            cout << "DEBUG: body_idxA = " << body_idxA << endl;
//            cout << "DEBUG: body_idxB = " << body_idxB << endl;

            // velocities and positions
            VectorXd u_A(6); VectorXd q_A(7); VectorXd u_B(6); VectorXd q_B(7);
            if (body_idxA != -1)
            {
                u_A = u_n.block<6,1>(6*body_idxA,0);
                q_A = q_n.block<7,1>(7*body_idxA,0);
            }
            else
            {
                u_A.setZero();
                q_A << 0, 0 , 0 , 1 , 0 ,0, 0;
            }
            if (body_idxB != -1)
            {
                u_B = u_n.block<6,1>(6*body_idxB,0);
                q_B = q_n.block<7,1>(7*body_idxB,0);
            }
            else {
                u_B.setZero();
                q_B << 0, 0, 0, 1, 0, 0, 0;
            }

            // collision values
            Vector3d normal = contact->n;
            Vector3d contact_wrt_bodyA = contact->contact_wrt_bodyA;
            Vector3d contact_wrt_bodyB = contact->contact_wrt_bodyB;


            // body A configuration
            Vector3d pos_A = q_A.head(3);
            Quaterniond rot_A(q_A(3),q_A(4),q_A(5),q_A(6));
            rot_A.normalize();

            // body B configuration
            Vector3d pos_B = q_B.head(3);
            Quaterniond rot_B(q_B(3),q_B(4),q_B(5),q_B(6));
            rot_B.normalize();

            // get the contact pts
            Vector3d contact_ptA = rot_A * contact_wrt_bodyA + pos_A;
//            cout << "DEBUG: pos_A = " << pos_A.transpose() << endl;
//            cout << "DEBUG: contact_wrt_bodyA = " << contact_wrt_bodyA.transpose() << endl;
//            cout << "DEBUG: contact_ptA = " << contact_ptA.transpose() << endl;
            Vector3d contact_ptB = rot_B * contact_wrt_bodyB + pos_B;
//            cout << "DEBUG: pos_B = " << pos_B.transpose() << endl;
//            cout << "DEBUG: contact_wrt_bodyB = " << contact_wrt_bodyB.transpose() << endl;
//            cout << "DEBUG: contact_ptB = " << contact_ptB.transpose() << endl;

            // compute penetration depth (should alread be done but w/e)
            double penetration_depth = (contact_ptA - contact_ptB).dot(normal);

            // if the penetration depth is > 0 constraint is inactive
            if (penetration_depth > 0) {
                is_active[cidx] = 0;
            } else {
                is_active[cidx] = 1;
            }

            // the constraint value
            double C_n = penetration_depth;

            // the tangent basis
            MatrixXd D = tangent_basis(normal);

            if (body_idxA != -1)
            {
                // the kinematic map matrices
                MatrixXd G_A = RigidObject::kinematic_map_G(rot_A);
//                cout << "DEBUG: G_A = " << G_A << endl;

                // gradients of the constraint for normal
                VectorXd grad_C_n_A(6); grad_C_n_A.head(3) = normal; grad_C_n_A.tail(3) = vec2skew(-normal).transpose()*contact_wrt_bodyA;
//                cout << "DEBUG: grad_C_n_A  = " << grad_C_n_A.transpose() << endl;

                ///////////////////////////////////////////////////////////////////////
                ///////// Assemble the Normal Constraint Matrices for Object A ///////
                //////////////////////////////////////////////////////////////////////

                VectorXd phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(cidx)) * grad_C_n_A;
//                cout << "DEBUG: phi_normal_wrt_q = " << phi_normal_wrt_q << endl;
                assert(!phi_normal_wrt_q.hasNaN());
                double phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(cidx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(cidx, cidx) = phi_normal_wrt_lambda;
                J_normal.block<1, 6>(cidx, 6*body_idxA) = phi_normal_wrt_q.transpose();

                // fill the normal constraint vector and normal jacobian
                phi_normal(cidx) = fischer(C_n, r_normal_ * lambda_normal(cidx));

                // Friction basis
                MatrixXd D_n_A = D;
                D_n_A.conservativeResize(6,2);
                D_n_A.block<3,1>(3,0)= vec2skew(-D_n_A.block<3,1>(0,0)).transpose()*(contact_wrt_bodyA);
                D_n_A.block<3,1>(3,1)= vec2skew(-D_n_A.block<3,1>(0,1)).transpose()*(contact_wrt_bodyA);

                ////////////////////////////////////////////////////////////////////////
                ///////// Assemble the Friction Constraint Matrices for Object A ///////
                ////////////////////////////////////////////////////////////////////////

                auto phi_friction_wrt_qdot = D_n_A;
                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
                double W_f = W_fischer(D_n_A, u_A, lambda_f, r_friction_, 0.1, lambda_normal(cidx));
                MatrixXd phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();

                // fill the friction vector
                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_A.transpose() * u_A + W_f * lambda_f;
            }

            if (body_idxB != -1)
            {
                // the kinematic map matrices
                MatrixXd G_B = RigidObject::kinematic_map_G(rot_B);
//                cout << "DEBUG: G_B = " << G_B << endl;

                // gradients of the constraint for normal
                VectorXd grad_C_n_B(6); grad_C_n_B.head(3) = -normal; grad_C_n_B.tail(3) = vec2skew(normal).transpose()*contact_wrt_bodyB;
//                cout << "DEBUG: grad_C_n_B  = " << grad_C_n_B.transpose() << endl;

                ///////////////////////////////////////////////////////////////////////
                ///////// Assemble the Normal Constraint Matrices for Object B ///////
                //////////////////////////////////////////////////////////////////////

                VectorXd phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(cidx)) * grad_C_n_B;
                assert(!phi_normal_wrt_q.hasNaN());
                double phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(cidx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(cidx, cidx) = phi_normal_wrt_lambda;
                J_normal.block<1, 6>(cidx, 6*body_idxB) = phi_normal_wrt_q.transpose();

                // fill the normal constraint vector and normal jacobian
                phi_normal(cidx) = fischer(C_n, r_normal_ * lambda_normal(cidx));

                // Friction basis
                MatrixXd D_n_B = -D;
                D_n_B.conservativeResize(6,2);
                D_n_B.block<3,1>(3,0)= vec2skew(-D_n_B.block<3,1>(0,0)).transpose()*contact_wrt_bodyB;
                D_n_B.block<3,1>(3,1)= vec2skew(-D_n_B.block<3,1>(0,1)).transpose()*contact_wrt_bodyB;

                ///////////////////////////////////////////////////////////////////////
                ///////// Assemble the Friction Constraint Matrices for Object B ///////
                //////////////////////////////////////////////////////////////////////

                auto phi_friction_wrt_qdot = D_n_B;
                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
                double W_f = W_fischer(D_n_B, u_B, lambda_f, r_friction_, 0.1, lambda_normal(cidx));
                MatrixXd phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();

                // fill the friction vector
                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_B.transpose() * u_B + W_f * lambda_f;
            }

        }

        cout << "DEBUG: is_active = " << is_active.transpose() << endl;

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = M_tilde;
//        cout << "DEBUG : H = " << H << endl;
        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());
        MatrixXd J(J_normal.rows() + J_friction.rows(), J_normal.cols());
        J.setZero();
        J << J_normal, J_friction;
        assert(!J.hasNaN());

        MatrixXd C(S.rows() + W.rows(),S.cols() + W.cols());
        C.setZero();
        C.block(0, 0, S.rows(), S.cols()) = S / pow(dt, 2);
        C.block(S.rows(), S.cols(), W.rows(), W.cols()) = W / dt;
        assert(!C.hasNaN());

        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(M_tilde.cols() == u_n.size());
        VectorXd g = M_tilde * (u_n) - dt * J.transpose() * lambda;
//        cout << "DEBUG: u " << u_n.transpose() << " , u_tilde = " << u_tilde.transpose() << endl;
        assert(!g.hasNaN());
//        cout << "DEBUG : g = " << g.transpose() << endl;
        MatrixXd I_n(phi_normal.size(), phi_normal.size());
        I_n.setIdentity();
        I_n /= dt;
        assert(!I_n.hasNaN());
        MatrixXd I_f(phi_friction.size(), phi_friction.size());
        I_f.setIdentity();
        assert(!I_f.hasNaN());
        MatrixXd I(I_n.rows() + I_f.rows(), I_n.cols() + I_f.cols());
        I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
        I.block(I_n.rows(), I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
        VectorXd phi(phi_normal.size() + phi_friction.size());
        phi << phi_normal, phi_friction;
        assert(!phi.hasNaN());
        VectorXd h = I * phi;

        assert(!h.hasNaN());
        assert(!J.hasNaN());
        assert(J.rows() == C.rows());
        assert(J.rows() == C.cols());
        auto H_inv = H.inverse();
        assert(!H_inv.hasNaN());
        MatrixXd A = J * H_inv * J.transpose() + C;

        if (is_ill_conditioned(A,1000))
        {
            cout << "DEBUG: Matrix A is ill conditioned!" << endl;
        }

        assert(!A.hasNaN());
        assert(J.cols() == H.rows() && H.cols() == g.size());
        assert(J.rows() == h.size());
        VectorXd b = (J * H_inv * g - h) / dt;
        assert(!b.hasNaN());

        for (int contact = 0; contact < numContacts; contact++)
        {
            if (is_active[contact] == 0 )
            {
                // zero the normal
                A.row(contact) = RowVectorXd::Zero(A.cols());
                b(contact) = 0;

                // zero the friction
                A.row(numContacts + 2*contact) = RowVectorXd::Zero(A.cols());
                A.row(numContacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
                b(numContacts + 2*contact) = 0;
                b(numContacts + 2*contact + 1) = 0;
            }
        }

        // solve the system
        VectorXd delta_lambda = A.fullPivLu().solve(b);

        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;

        // solve this system with Gauss-Seidel
        assert(!delta_lambda.hasNaN());
        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
        cout << "DEBUG: delta_u = " << delta_u.transpose() << endl;
        assert(!delta_u.hasNaN());
        // perform a line search
        double t = 0.1;

        // update solution
        lambda += t*delta_lambda;
        u_n += t*delta_u;
        q_n += G*u_n;

        if (delta_lambda.squaredNorm() < 1e-12)
        {
            cout << "Newtons method converged in " << n << "iterations" << endl;
            return;
        }

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
        lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);

        r_normal = pow(dt,2)*(J*M_tilde_inv*J.transpose()).diagonal();
        r_friction = dt*(J*M_tilde_inv*J.transpose()).diagonal();
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////


void ImplicitLCP::setup_implicit_lcp_deformable_linear(RigidObject* rigidObject, DeformableObject* deformableObject, vector<Contact*> contacts,
                                                       double dt, int newton_it, int linear_it, double newton_tol, double linear_tol)
{

//    int ndof = 6 + 3*deformableObject->m_numVerts;
//    int numContacts = contacts.size();
//    int numTets = deformableObject->tetrahedra().rows();
//    VectorXd phi_elastic(3*numTets); phi_elastic.setZero();
//    VectorXd phi_normal(numContacts); phi_normal.setZero();
//    VectorXd phi_friction(2*numContacts); phi_friction.setZero();
//    VectorXd lambda_elastic(6*numTets); lambda_elastic.setZero();
//    VectorXd lambda_normal(numContacts); lambda_normal.setZero();
//    VectorXd lambda_friction(2*numContacts); lambda_friction.setZero();
//    VectorXd lambda(lambda_elastic.size() + lambda_normal.size() + lambda_friction.size()); lambda.setZero();
//    MatrixXd J_elastic(numTets,3); J_elastic.setZero();
//    MatrixXd J_normal(numContacts,3); J_normal.setZero();
//    MatrixXd J_friction(2*numContacts, 3); J_friction.setZero();
//    MatrixXd S(numContacts,numContacts); S.setZero();
//    MatrixXd W(2*numContacts,2*numContacts); W.setZero();
//    VectorXi is_active(numContacts);
//    VectorXd r_normal = VectorXd::Ones(numContacts);
//    VectorXd r_friction = VectorXd::Ones(numContacts);
//
//    for (int n = 0 ; n < newton_it; n++) {
//
//        // create the mass matrix
//        MatrixXd M(ndof,ndof);
//        M.block<6,6>(0,0) = rigidObject->mass_matrix();
//        M.block(6,6,3*deformableObject->m_numVerts,3*deformableObject->m_numVerts) =
//                deformableObject->mass_matrix();
//        MatrixXd M_inv(ndof,ndof);
//        M_inv.block<6,6>(0,0) = rigidObject->inverse_mass_matrix();
//        M_inv.block(6,6,3*deformableObject->m_numVerts,3*deformableObject->m_numVerts) =
//                deformableObject->inverse_mass_matrix();
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        //////////////////////////////////////// SETUP CONTACT CONSTRAINTS ////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//
//        // assemble the matrix
//        for (int cidx = 0; cidx < numContacts; cidx++) {
//
//            // set the preconditioner values
//            double r_friction_ = dt; double r_normal_ = dt*dt;
//
//            // the collision
//            auto &contact = contacts[cidx];
//
//            // set the contact index
//            int body_idxA = contact->bodyIdxA;
//            int body_idxB = contact->bodyIdxB;
//
//            // velocities and positions
//            VectorXd u_A(6); VectorXd q_A(7); VectorXd u_B(6); VectorXd q_B(7);
//            if (body_idxA != -1)
//            {
//                u_A = u_n.block<6,1>(6*body_idxA,0);
//                q_A = q_n.block<7,1>(7*body_idxA,0);
//            }
//            else
//            {
//                u_A.setZero();
//                q_A << 0, 0 , 0 , 1 , 0 ,0, 0;
//            }
//            if (body_idxB != -1)
//            {
//                u_B = u_n.block<6,1>(6*body_idxB,0);
//                q_B = q_n.block<7,1>(7*body_idxB,0);
//            }
//            else {
//                u_B.setZero();
//                q_B << 0, 0, 0, 1, 0, 0, 0;
//            }
//
//            // collision values
//            Vector3d normal = contact->normal;
//            Vector3d contact_wrt_bodyA = contact->contact_wrt_bodyA;
//            Vector3d contact_wrt_bodyB = contact->contact_wrt_bodyB;
//
//
//            // body A configuration
//            Vector3d pos_A = q_A.head(3);
//            Quaterniond rot_A(q_A(3),q_A(4),q_A(5),q_A(6));
//            rot_A.normalize();
//
//            // body B configuration
//            Vector3d pos_B = q_B.head(3);
//            Quaterniond rot_B(q_B(3),q_B(4),q_B(5),q_B(6));
//            rot_B.normalize();
//
//            // get the contact pts
//            Vector3d contact_ptA = rot_A * contact_wrt_bodyA + pos_A;
//            Vector3d contact_ptB = rot_B * contact_wrt_bodyB + pos_B;
//
//            // compute penetration depth (should alread be done but w/e)
//            double penetration_depth = (contact_ptA - contact_ptB).dot(normal);
//
//            // if the penetration depth is > 0 constraint is inactive
//            if (penetration_depth > 0) {
//                is_active[cidx] = 0;
//            } else {
//                is_active[cidx] = 1;
//            }
//
//            // the constraint value
//            double C_n = penetration_depth;
//
//            // the tangent basis
//            MatrixXd D = tangent_basis(normal);
//
//            if (body_idxA != -1)
//            {
//                // the kinematic map matrices
//                MatrixXd G_A = RigidObject::kinematic_map_G(rot_A);
//
//                // gradients of the constraint for normal
//                VectorXd grad_C_n_A(6); grad_C_n_A.head(3) = normal; grad_C_n_A.tail(3) = vec2skew(-normal).transpose()*contact_wrt_bodyA;
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Normal Constraint Matrices for Object A ///////
//                //////////////////////////////////////////////////////////////////////
//
//                VectorXd phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(cidx)) * grad_C_n_A;
////                cout << "DEBUG: phi_normal_wrt_q = " << phi_normal_wrt_q << endl;
//                assert(!phi_normal_wrt_q.hasNaN());
//                double phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(cidx)) * r_normal_;
//                assert(!isnan(phi_normal_wrt_lambda));
//                S(cidx, cidx) = phi_normal_wrt_lambda;
//                J_normal.block<1, 6>(cidx, 6*body_idxA) = phi_normal_wrt_q.transpose();
//
//                // fill the normal constraint vector and normal jacobian
//                phi_normal(cidx) = fischer(C_n, r_normal_ * lambda_normal(cidx));
//
//                // Friction basis
//                MatrixXd D_n_A = D;
//                D_n_A.conservativeResize(6,2);
//                D_n_A.block<3,1>(3,0)= vec2skew(-D_n_A.block<3,1>(0,0)).transpose()*(contact_wrt_bodyA);
//                D_n_A.block<3,1>(3,1)= vec2skew(-D_n_A.block<3,1>(0,1)).transpose()*(contact_wrt_bodyA);
//
//                ////////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Friction Constraint Matrices for Object A ///////
//                ////////////////////////////////////////////////////////////////////////
//
//                auto phi_friction_wrt_qdot = D_n_A;
//                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
//                double W_f = W_fischer(D_n_A, u_A, lambda_f, r_friction_, 0.1, lambda_normal(cidx));
//                MatrixXd phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();
//
//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_A.transpose() * u_A + W_f * lambda_f;
//            }
//
//            if (body_idxB != -1)
//            {
//                // the kinematic map matrices
//                MatrixXd G_B = RigidObject::kinematic_map_G(rot_B);
//
//                // gradients of the constraint for normal
//                VectorXd grad_C_n_B(6); grad_C_n_B.head(3) = -normal; grad_C_n_B.tail(3) = vec2skew(normal).transpose()*contact_wrt_bodyB;
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Normal Constraint Matrices for Object B ///////
//                //////////////////////////////////////////////////////////////////////
//
//                VectorXd phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(cidx)) * grad_C_n_B;
//                assert(!phi_normal_wrt_q.hasNaN());
//                double phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(cidx)) * r_normal_;
//                assert(!isnan(phi_normal_wrt_lambda));
//                S(cidx, cidx) = phi_normal_wrt_lambda;
//                J_normal.block<1, 6>(cidx, 6*body_idxB) = phi_normal_wrt_q.transpose();
//
//                // fill the normal constraint vector and normal jacobian
//                phi_normal(cidx) = fischer(C_n, r_normal_ * lambda_normal(cidx));
//
//                // Friction basis
//                MatrixXd D_n_B = -D;
//                D_n_B.conservativeResize(6,2);
//                D_n_B.block<3,1>(3,0)= vec2skew(-D_n_B.block<3,1>(0,0)).transpose()*contact_wrt_bodyB;
//                D_n_B.block<3,1>(3,1)= vec2skew(-D_n_B.block<3,1>(0,1)).transpose()*contact_wrt_bodyB;
//
//                ///////////////////////////////////////////////////////////////////////
//                ///////// Assemble the Friction Constraint Matrices for Object B ///////
//                //////////////////////////////////////////////////////////////////////
//
//                auto phi_friction_wrt_qdot = D_n_B;
//                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * cidx, 0);
//                double W_f = W_fischer(D_n_B, u_B, lambda_f, r_friction_, 0.1, lambda_normal(cidx));
//                MatrixXd phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * cidx, 2 * cidx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 6>(2 * cidx, 0) = phi_friction_wrt_qdot.transpose();
//
//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * cidx, 0) = D_n_B.transpose() * u_B + W_f * lambda_f;
//            }
//
//        }
//
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        //////////////////////////////////////// SETUP ELASTIC CONSTRAINTS ////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
//        for (int tidx = 0; tidx < numTets; tidx++)
//        {
//            // FOR EACH TETRAHEDRA COMPUTE THE CONSTRAINT FORCE
//            int v0 = tets(tidx,0); int v1 = tets(tidx,1); int v2 = tets(tidx,2); int v3 = tets(tidx,3);
//        }
//
//        assert(!J_friction.hasNaN());
//        assert(!J_normal.hasNaN());
//        assert(!S.hasNaN());
//        assert(!W.hasNaN());
//        MatrixXd H = M_tilde;
//        assert(!H.hasNaN());
//        assert(J_normal.cols() == J_friction.cols());
//        MatrixXd J(J_normal.rows() + J_friction.rows(), J_normal.cols());
//        J.setZero();
//        J << J_normal, J_friction;
//        assert(!J.hasNaN());
//
//        MatrixXd C(S.rows() + W.rows(),S.cols() + W.cols());
//        C.setZero();
//        C.block(0, 0, S.rows(), S.cols()) = S / pow(dt, 2);
//        C.block(S.rows(), S.cols(), W.rows(), W.cols()) = W / dt;
//        assert(!C.hasNaN());
//
//        assert(!lambda.hasNaN());
//        assert(J.rows() == lambda.size());
//        assert(M_tilde.cols() == u_n.size());
//        VectorXd g = M_tilde * (u_n) - dt * J.transpose() * lambda;
//        assert(!g.hasNaN());
//        MatrixXd I_b = MatrixXd::Identity(6*numTets,6*numTets);
//        I_b /= dt;
//        MatrixXd I_n = MatrixXd::Identity(phi_normal.size(),phi_normal.size());
//        I_n /= dt;
//        assert(!I_n.hasNaN());
//        MatrixXd I_f = MatrixXd::Identity(phi_friction.size(), phi_friction.size());
//        assert(!I_f.hasNaN());
//        MatrixXd I(I_n.rows() + I_f.rows(), I_n.cols() + I_f.cols());
//        I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
//        I.block(I_n.rows(), I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
//        VectorXd phi(phi_normal.size() + phi_friction.size());
//        phi << phi_normal, phi_friction;
//        assert(!phi.hasNaN());
//        VectorXd h = I * phi;
//
//        assert(!h.hasNaN());
//        assert(!J.hasNaN());
//        assert(J.rows() == C.rows());
//        assert(J.rows() == C.cols());
//        auto H_inv = H.inverse();
//        assert(!H_inv.hasNaN());
//        MatrixXd A = J * H_inv * J.transpose() + C;
//
//        if (is_ill_conditioned(A,1000))
//        {
//            cout << "DEBUG: Matrix A is ill conditioned!" << endl;
//        }
//
//        assert(!A.hasNaN());
//        assert(J.cols() == H.rows() && H.cols() == g.size());
//        assert(J.rows() == h.size());
//        VectorXd b = (J * H_inv * g - h) / dt;
//        assert(!b.hasNaN());
//
//        for (int contact = 0; contact < numContacts; contact++)
//        {
//            if (is_active[contact] == 0 )
//            {
//                // zero the normal
//                A.row(contact) = RowVectorXd::Zero(A.cols());
//                b(contact) = 0;
//
//                // zero the friction
//                A.row(numContacts + 2*contact) = RowVectorXd::Zero(A.cols());
//                A.row(numContacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
//                b(numContacts + 2*contact) = 0;
//                b(numContacts + 2*contact + 1) = 0;
//            }
//        }
//
//        // solve the system
//        VectorXd delta_lambda = A.fullPivLu().solve(b);
//
//        // solve this system with Gauss-Seidel
//        assert(!delta_lambda.hasNaN());
//        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
//        assert(!delta_u.hasNaN());
//        // perform a line search
//        double t = 0.1;
//
//        // update solution
//        lambda += t*delta_lambda;
//        u_n += t*delta_u;
//        q_n += G*u_n;
//
//        if (delta_lambda.squaredNorm() < 1e-12)
//        {
//            cout << "Newtons method converged in " << n << "iterations" << endl;
//            return;
//        }
//
//        // set the constraint forces
//        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
//        lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);
//
//        r_normal = pow(dt,2)*(J*M_tilde_inv*J.transpose()).diagonal();
//        r_friction = dt*(J*M_tilde_inv*J.transpose()).diagonal();
//    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////

// whats the point in doing a quasistatic approximation
void ImplicitLCP::setup_implicit_lcp_haptic_quasistatic(VectorXd& q_n, const Eigen::VectorXd &q_tilde, //! decide whether this should update (device)
                                                        Eigen::VectorXd &u_n, //! also determine which velocity to send
                                                        vector<Contact*> &contacts, double dt, int newton_it, int linear_it,
                                                        double newton_tol, double linear_tol)
{
    assert(!q_tilde.hasNaN()); assert(!q_n.hasNaN());

    // linear translation
    int ndof = 3;

    MatrixXd K = Matrix3d().setIdentity();
    MatrixXd K_inv = K.inverse();
    int numContacts = contacts.size();
    VectorXd phi_normal(numContacts); phi_normal.setZero();
    VectorXd phi_friction(2*numContacts); phi_friction.setZero();
    VectorXd lambda_normal(numContacts); lambda_normal.setZero();
    VectorXd lambda_friction(2*numContacts); lambda_friction.setZero();
    VectorXd lambda(lambda_normal.size()); lambda.setZero();
//    VectorXd lambda(lambda_normal.size() + lambda_friction.size()); lambda.setZero();
    MatrixXd J_normal(numContacts,ndof); J_normal.setZero();
    MatrixXd J_friction(2*numContacts,ndof); J_friction.setZero();
    MatrixXd E;
    MatrixXd S(numContacts,numContacts); S.setZero();
    MatrixXd W(2*numContacts,2*numContacts); W.setZero();
    VectorXi is_active(numContacts);
    VectorXd r_normal = VectorXd::Ones(numContacts);
    VectorXd r_friction = VectorXd::Ones(numContacts);

    for (int n = 0 ; n < newton_it; n++) {

        cout << "Newton It : " << n << endl;
        // assemble the matrix
        for (int colIdx = 0; colIdx < numContacts; colIdx++) {

            // set the preconditioner values
            double r_friction_ = r_friction(colIdx); double r_normal_ = r_normal(colIdx);

            // collision values
            Vector3d normal = contacts[colIdx]->n;
            Vector3d contact_wrt_bodyA = contacts[colIdx]->contact_wrt_bodyA;
            Vector3d contact_wrt_bodyB = contacts[colIdx]->contact_wrt_bodyB;

            // set the contact index
            int contact_idx_bodyA = contacts[colIdx]->bodyIdxA;
            int contact_idx_bodyB = contacts[colIdx]->bodyIdxB;

            // get the contact pts
            //! need to do the full rotation transform
            Vector3d contact_ptA;
            if (contact_idx_bodyA != -1)
            {
                contact_ptA = q_n.block<3, 1>(contact_idx_bodyA, 0) + contact_wrt_bodyA;
            }
            else
            {
                contact_ptA =  contact_wrt_bodyA;
            }
            Vector3d contact_ptB;
            if (contact_idx_bodyB != -1 )
            {
                contact_ptB = q_n.block<3,1>(contact_idx_bodyB,0) + contact_wrt_bodyB;
            }
            else{
                contact_ptB =  contact_wrt_bodyB;
            }


            double penetration_depth = (contact_ptA - contact_ptB).dot(normal);

            if (penetration_depth > 0 )
            {
                is_active[colIdx] = 0;
            }
            else
            {
                is_active[colIdx] = 1;
            }

            double C_n = penetration_depth;

            VectorXd grad_C_n_A = normal.normalized();
            VectorXd grad_C_n_B = -normal.normalized();

            //
            if (contact_idx_bodyA == 0)
            {
                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_A;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = checkNearZero(phi_normal_wrt_q).transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));

                // now the same for the friction jacobian
//                auto D = tangent_basis(normal);
//                auto phi_friction_wrt_qdot = D;
//                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * colIdx, 0);
//                double W_f = W_fischer(D, u_n, lambda_f, r_friction_, 0.1, lambda_normal(colIdx));
//                auto phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * colIdx, 2 * colIdx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 3>(2 * colIdx, 0) = phi_friction_wrt_qdot.transpose() * G;

//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * colIdx, 0) = D.transpose() * u_n + W_f * lambda_f;
            }
            if (contact_idx_bodyB == 0)
            {
                auto phi_normal_wrt_q = fischer_partial_a(C_n, r_normal_ * lambda_normal(colIdx)) * grad_C_n_B;
                assert(!phi_normal_wrt_q.hasNaN());
                auto phi_normal_wrt_lambda = fischer_partial_b(C_n, r_normal_ * lambda_normal(colIdx)) * r_normal_;
                assert(!isnan(phi_normal_wrt_lambda));
                S(colIdx, colIdx) = phi_normal_wrt_lambda;
                auto G = MatrixXd(3, 3).setIdentity(); // Replace later with transform
                J_normal.block<1, 3>(colIdx, 0) = checkNearZero(phi_normal_wrt_q).transpose() * G;

                // fill the normal constraint vector and normal jacobian
                phi_normal(colIdx) = fischer(C_n, r_normal_ * lambda_normal(colIdx));

                // now the same for the friction jacobian
//                auto D = tangent_basis(normal);
//                auto phi_friction_wrt_qdot = D;
//                VectorXd lambda_f = lambda_friction.block<2, 1>(2 * colIdx, 0);
//                double W_f = W_fischer(D, u_n, lambda_f, r_friction_, 0.1, lambda_normal(colIdx));
//                auto phi_friction_wrt_lambda_f = W_f * Matrix2d::Identity();
//                W.block<2, 2>(2 * colIdx, 2 * colIdx) = phi_friction_wrt_lambda_f;
//                J_friction.block<2, 3>(2 * colIdx, 0) = phi_friction_wrt_qdot.transpose() * G;
//
//                // fill the friction vector
//                phi_friction.block<2, 1>(2 * colIdx, 0) = D.transpose() * u_n + W_f * lambda_f;
            }
        }

//        cout << "DEBUG: S = " << S << endl;
//        cout << "DEBUG : phi_normal = " << phi_normal.transpose() << endl;
//        cout << "DEBUG : J_normal = " << J_normal << endl;

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = K;
//        cout << "DEBUG : H = " << H << endl;
        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());
        MatrixXd J(J_normal.rows(), J_normal.cols());
//        J.resize(J_normal.rows() + J_friction.rows(), J_normal.cols());
        J.setZero();
        J << J_normal;
//        J << J_normal, J_friction;
        assert(!J.hasNaN());

        MatrixXd C(E.rows() + S.rows(), E.cols() + S.cols());
//        C.resize(E.rows() + S.rows() + W.rows(), E.cols() + S.cols() + W.cols());
        C.setZero();
        C.block(0, 0, E.rows(), E.cols()) = E / pow(dt, 2);
        C.block(E.rows(), E.cols(), S.rows(), S.cols()) = S / dt;
//        C.block(E.rows() + S.rows(), E.cols() + S.cols(), W.rows(), W.cols()) = W / dt;
        assert(!C.hasNaN());

        VectorXd delta_lambda(lambda_normal.size());
//        VectorXd delta_lambda(lambda_normal.size() + lambda_friction.size());
        delta_lambda.setZero();
        // not sure if correct
        // delta_lambda << lambda_normal, lambda_friction;
        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(K.cols() == u_n.size());
        VectorXd g = J.transpose() * lambda;
        assert(!g.hasNaN());
        MatrixXd I_n(phi_normal.size(), phi_normal.size());
        I_n.setIdentity();
        I_n /= dt;
        assert(!I_n.hasNaN());
        MatrixXd I_f(phi_friction.size(), phi_friction.size());
        I_f.setIdentity();
        assert(!I_f.hasNaN());
        MatrixXd I(I_n.rows(), I_n.cols());
//        MatrixXd I(I_n.rows() + I_f.rows(), I_n.cols() + I_f.cols());
        I.block(0, 0, I_n.rows(), I_n.cols()) = I_n;
//        cout << "DEBUG : I = " << I << endl;
//        I.block(I_n.rows(), I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
        VectorXd phi(phi_normal.size());
//        VectorXd phi(phi_normal.size() + phi_friction.size());
        phi << phi_normal;
//        phi << phi_normal, phi_friction;
        assert(!phi.hasNaN());
//        cout << "DEBUG: phi = " << phi.transpose() << endl;
        VectorXd h = I * phi;
//        cout << "DEBUG: h = " << h.transpose() << endl;

        assert(!h.hasNaN());
        assert(!J.hasNaN());
        assert(J.rows() == C.rows());
        assert(J.rows() == C.cols());
        auto H_inv = H.inverse();
        assert(!H_inv.hasNaN());
        MatrixXd A = J * H_inv * J.transpose() + C;

        if (is_ill_conditioned(A,1000))
        {
            cout << "DEBUG: Matrix A is ill conditioned!" << endl;
        }

        assert(!A.hasNaN());
        assert(J.cols() == H.rows() && H.cols() == g.size());
        assert(J.rows() == h.size());
        VectorXd b = (J * H_inv * g - h);
//        cout << "DEBUG: b = " << b.transpose() << endl;
//        cout << "DEBUG: A = " << A << endl;
        assert(!b.hasNaN());

        for (int contact = 0; contact < numContacts; contact++)
        {
            if (is_active[contact] == 0 )
            {
                // zero the normal
//                A.row(contact) = RowVectorXd::Zero(A.cols());
//                b(contact) = 0;

//                A.row(numContacts + 2*contact) = RowVectorXd::Zero(A.cols());
//                A.row(numContacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
//                b(numContacts + 2*contact) = 0;
//                b(numContacts + 2*contact + 1) = 0;
            }
        }

        cout << "DEBUG: is_active = " << is_active.transpose() << endl;
        // solve this system with Gauss-Seidel
        auto res = gauss_seidel(A, b, delta_lambda,linear_tol);
        cout << "DEBUG: delta_lambda = " << delta_lambda.transpose() << endl;
        assert(!delta_lambda.hasNaN());
        auto delta_q = H_inv*(J.transpose()*delta_lambda - g);
        cout << "DEBUG: delta_q = " << delta_q.transpose() << endl;
        assert(!delta_q.hasNaN());
        // perform a line search
        double t = 0.75;
        // update solution
        lambda += t*delta_lambda;
        u_n += t*delta_q/dt;
        q_n += t*delta_q;

        cout << "DEBUG: q_n = " << q_n.transpose() << endl;

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
//            lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);

        auto f_normal = J_normal.transpose()*lambda_normal;
        auto f_friction = J_friction.transpose()*lambda_friction;

        // if the residual is low enough then break out
        r_normal = pow(dt,2)*(J*K_inv*J.transpose()).diagonal();
        r_friction = dt*(J*K_inv*J.transpose()).diagonal();
    }

//    cout << "Newtons Method Failed to Converge" << endl;
}