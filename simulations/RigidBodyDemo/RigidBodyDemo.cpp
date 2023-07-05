//
// Created by root on 7/4/23.
//

#include "RigidBodyDemo.h"


void RigidBodyDemo::initialize()
{
    // add balls to the scene
    for (int i = 0 ; i < 4; i++)
    {
        if (i == 0)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj"));
            m_objects.back()->set_local_pos(Vector3d(0,sqrt(3)/4*0.55,0.0));
            m_objects.back()->getMesh(0)->m_material->setBlue();
        }
        else if(i == 1)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj"));
            m_objects.back()->set_local_pos(Vector3d(-0.55/2,-sqrt(3)/4*0.55,0.0));
            m_objects.back()->getMesh(0)->m_material->setRed();
        }
        else if (i == 2)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj"));
            m_objects.back()->set_local_pos(Vector3d(0.55/2,-sqrt(3)/4*0.55,0.0));
            m_objects.back()->getMesh(0)->m_material->setYellow();
        }
        else if (i == 3)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/RigidBodyDemo/bowl.obj"));
            m_objects.back()->set_is_static(true);
            m_objects.back()->getMesh(0)->m_material->setBrownCornsilk();
            m_objects.back()->set_local_rot(generateQuaternion(90,0,0));
            m_objects.back()->scale(7.5);
        }
        m_world->addChild(m_objects.back());
        m_objects.back()->import_mesh_data();
    }

    // create the collision detector
    m_collisionDetector = new CollisionDetector(m_objects);
}

void RigidBodyDemo::step(double dt)
{

    // reset forces
    for(auto b : m_objects)
    {
        b->f = b->mass * Vector3d(0., 0, -0.0098);
        b->tau.setZero();
        b->fc.setZero();
        b->tauc.setZero();
        b->m_contacts.clear();
    }

    // computes any collisions
    m_collisionDetector->computeCollisions();

    // computes the constraints
    if (m_collisionDetector->m_contacts.size() > 0)
    {
        computeConstraints(dt);
    }

    // explicit step
    for(auto* b : m_objects)
    {
        b->x_last = b->x; b->q_last = b->q_last; b->xdot_last = b->xdot; b->omega_last = b->omega;
        if( !b->is_static )
        {
            b->xdot += dt * (1./b->mass) * (b->f + b->fc);
            b->omega += dt * b->Iinv * (b->tau + b->tauc);
            b->x += dt * b->xdot;
            b->q = sumQuaternions(b->q,
                                        multiplyQuaternionScalar(dt * 0.5 , Quaterniond(0, b->omega[0], b->omega[1], b->omega[2]) * b->q));
            b->q.normalize();
        }
        else
        {
            b->xdot.setZero();
            b->omega.setZero();
        }
    }
}

void RigidBodyDemo::updateGraphics()
{
    for (auto object : m_objects)
    {
        object->update_mesh_position();
    }
}

void RigidBodyDemo::computeConstraints(double dt)
{
    // get contacts
    auto contacts = m_collisionDetector->m_contacts;
    // get the number of objects
    int num_objects = m_objects.size();
    // get the number of contatcs
    int num_contacts = contacts.size();

    VectorXd phi_normal(num_contacts); phi_normal.setZero();
    VectorXd phi_friction(2 * num_contacts); phi_friction.setZero();
    VectorXd lambda_normal(num_contacts); lambda_normal.setZero();
    VectorXd lambda_friction(2 * num_contacts); lambda_friction.setZero();
    VectorXd lambda(lambda_normal.size() + lambda_friction.size()); lambda.setZero();
    MatrixXd J_normal(num_contacts, num_objects * 6); J_normal.setZero();
    MatrixXd J_friction(2 * num_contacts, num_objects * 6); J_friction.setZero();
    MatrixXd S(num_contacts, num_contacts); S.setZero();
    MatrixXd W(2 * num_contacts, 2 * num_contacts); W.setZero();
    VectorXi is_active(num_contacts); // determines if constraint is active or not

    VectorXd q(7 * num_objects); // the constrained position in global coordinates
    VectorXd u(6 * num_objects); // the constrained velocity in generalized coordinates
    VectorXd u_tilde(6 * num_objects); // the unconstrained velocity

    // fill the unconstrained positions
    for (int idx = 0 ; idx < num_objects; idx++)
    {
        u_tilde.block<3,1>(6*idx,0) = m_objects[idx]->xdot;
        u_tilde.block<3,1>(6*idx+3, 0) = m_objects[idx]->omega;
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////// START THE NEWTON ITERATIONS ////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int n = 0; n < m_settings.newtonIterations; n++)
    {
        MatrixXd M_tilde(6*num_objects,6*num_objects); M_tilde.setZero();
        MatrixXd M_inv_tilde(6*num_objects,6*num_objects); M_inv_tilde.setZero();
        for (int idx = 0; idx < num_objects; idx++)
        {
            m_objects[idx]->update_inertia_matrix();
            M_tilde.block<6,6>(6*idx,6*idx) = m_objects[idx]->generalized_mass();
            M_inv_tilde.block<6,6>(6*idx,6*idx) = m_objects[idx]->generalized_mass_inverse();
            u.block<6,1>(6*idx,0) = m_objects[idx]->generalized_vel();
        }

        for (int cidx = 0; cidx < num_contacts; cidx++) {

            // preconditioning values
            double r_n = dt*dt; double r_f = dt;

            // the collision
            auto &contact = contacts[cidx];

            // compute the constraint data
            contact->compute_constraints(r_n,r_f,
                                         lambda_normal(cidx),lambda_friction.block<2,1>(2*cidx,0));

            // object indices
            int objectA_idx = contact->objectA->m_idx;
            int objectB_idx = contact->objectB->m_idx;

            // contact data
            auto J_n_A = contact->Jacobian_normal_A(); auto J_n_B = contact->Jacobian_normal_B();
            auto J_f_A = contact->Jacobian_friction_A(); auto J_f_B = contact->Jacobian_friction_B();
            auto phi_n = contact->phi_normal();
            auto phi_f = contact->phi_friction();
            auto phi_n_wrt_lambda_n = contact->partial_phi_normal_wrt_lambda_normal();
            auto phi_f_wrt_lambda_f = contact->partial_phi_friction_wrt_lambda_friction();

            // fill the normal constraint vector and normal jacobian
            phi_normal(cidx) = phi_n;
            S(cidx, cidx) = phi_n_wrt_lambda_n;
            J_normal.block<1, 6>(cidx, 6*objectA_idx) = J_n_A;
            J_normal.block<1, 6>(cidx, 6*objectA_idx) = J_n_B;

            // fill the friction constraint vector and friction jacobian
            phi_friction.block<2, 1>(2 * cidx, 0) = phi_f;
            W.block<2, 2>(2 * cidx, 2 * cidx) = phi_f_wrt_lambda_f;
            J_friction.block<2, 6>(2 * cidx, 6*objectA_idx) = J_f_A;
            J_friction.block<2, 6>(2 * cidx, 6*objectB_idx) = J_f_B;

        }

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = M_tilde;
        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());

        MatrixXd J(J_normal.rows() + J_friction.rows(), J_normal.cols());
        J.setZero();
        J << J_normal, J_friction;
        assert(!J.hasNaN());

        MatrixXd C(S.rows() + W.rows(),S.cols() + W.cols());
        C.setZero();
        C.block(0, 0, S.rows(), S.cols()) = S / pow(dt, 2);
        C.block( S.rows(), S.cols(), W.rows(), W.cols()) = W / dt;
        assert(!C.hasNaN());
        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(M_tilde.cols() == u.size());
        VectorXd g = M_tilde * (u - u_tilde) - dt * J.transpose() * lambda;
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
        // zero out any nans caused by inversion
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

        // zero out any inactive contacts
        for (int contact = 0; contact < num_contacts; contact++)
        {
            if (is_active[contact] == 0 )
            {
                A.row(num_contacts + 2*contact) = RowVectorXd::Zero(A.cols());
                A.row(num_contacts + 2*contact + 1) = RowVectorXd::Zero(A.cols());
                b(num_contacts + 2*contact) = 0;
                b(num_contacts + 2*contact + 1) = 0;
            }
        }

        // solve this system with Gauss-Seidel
        VectorXd delta_lambda = A.fullPivLu().solve(b);

        assert(!delta_lambda.hasNaN());
        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
        assert(!delta_u.hasNaN());
        // perform a line search
        int t = 0.5;

        // update solution
        lambda += t*delta_lambda;
        u += t*delta_u;
        for (int bidx = 0 ; bidx < m_objects.size(); bidx++)
        {
            q.block<7,1>(bidx*7,0) += dt*m_objects[bidx]->kinematic_map_G()*u.block<7,1>(bidx*7,0);
            m_objects[bidx]->x = q.block<3,1>(bidx*7,0); m_objects[bidx]->q = q.block<4,1>(bidx*7+3,0);
            m_objects[bidx]->xdot = u.block<3,1>(bidx*6,0); m_objects[bidx]->omega = q.block<3,1>(bidx*6+3,0);
        }

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
        lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);


    }
}