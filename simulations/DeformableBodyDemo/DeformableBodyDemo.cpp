#include "DeformableBodyDemo.h"

void DeformableBodyDemo::initialize()
{
    // add a gummy bear to the scene
    int i = 0;
    m_objects.emplace_back(new DeformableObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/DeformableBodyDemo/beam.off"));
    m_objects.back()->set_local_pos(Vector3d(0,sqrt(3)/4*0.55,0.0));
    m_objects.back()->getMesh(0)->m_material->setPurple();
    m_objects.back()->getMesh(0)->m_material->setShininess(0.5);

    m_world->addChild(m_objects.back());
    m_objects.back()->import_mesh_data();

    // find the beam leftmost points and pin them
    auto& obj = m_objects.back();
    double xmin = obj->m_vertices.col(0).minCoeff();
    for (int vidx = 0 ; vidx < obj->m_numVerts; vidx++)
    {
        Vector3d p0 = obj->m_vertices.row(vidx);
        if (p0(0) < xmin + 1e-6) obj->fixed_constraints.emplace_back(new FixedConstraint(vidx,p0));
    }

    // create the linear elastic constraints
    auto& x = obj->x;
    for (int tidx = 0; tidx < obj->m_numTets; tidx++)
    {
        int v0 = obj->m_tetrahedra(tidx,0);
        int v1 = obj->m_tetrahedra(tidx,1);
        int v2 = obj->m_tetrahedra(tidx,2);
        int v3 = obj->m_tetrahedra(tidx,3);
        obj->linear_hookean_constraints.emplace_back(new LinearHookeanConstraint(v0,v1,v2,v3,tidx,x));
    }

    // get the size of the mass matrix
    int mdim = 0;
    for (auto obj : m_objects) mdim += 3*obj->m_numVerts;
    M.resize(mdim,mdim); M_inv.resize(mdim,mdim);
    int dimctr = 0;
    for (auto obj : m_objects)
    {
        M.block(dimctr,dimctr,3*obj->m_numVerts,3*obj->m_numVerts) = obj->mass_matrix();
        M_inv.block(dimctr,dimctr,3*obj->m_numVerts,3*obj->m_numVerts) = obj->inverse_mass_matrix();
        dimctr += 3*obj->m_numVerts;
    }

    // create the collision detector
    m_collisionDetector = new CollisionDetector(m_objects);

}

void DeformableBodyDemo::step(double dt)
{

    // reset forces
    for(auto b : m_objects)
    {
        for (int vidx = 0 ; vidx < b->m_numVerts; vidx++)
        b->f = b->mass * createGravityVector(b->m_numVerts);
        b->fc.setZero();
        b->m_contacts.clear();
    }

    // computes any collisions
    //m_collisionDetector->computeCollisions();

    // computes the constraints
    computeConstraints(dt);


    // explicit step
    for(auto* b : m_objects)
    {
        b->x_last = b->x;  b->xdot_last = b->xdot;

        if( !b->is_static )
        {
            b->xdot += dt * (1./b->mass) * (b->f + b->fc);
            b->x += dt * b->xdot;
        }
        else
        {
            b->xdot.setZero();
        }
    }
}

void DeformableBodyDemo::updateGraphics()
{
    for (auto object : m_objects)
    {
        object->update_mesh_position();
    }
}

void DeformableBodyDemo::computeConstraints(double dt)
{
    auto contacts = m_collisionDetector->m_contacts; // get contacts
    int num_objects = m_objects.size(); // get the number of objects
    int num_contacts = contacts.size(); // get the number of contatcs
    int num_bilateral_constraints = 0;
    int num_elements = 0;
    for (auto obj : m_objects)
    {
        num_bilateral_constraints += obj->m_numTets;
        num_elements += obj->m_numVerts;
    }

    VectorXd phi_normal(num_contacts); phi_normal.setZero(); //  the normal constraint values
    VectorXd phi_friction(2 * num_contacts); phi_friction.setZero(); //  the friction constraint values
    VectorXd h_bilateral(9 * num_bilateral_constraints); h_bilateral.setZero(); // the bilateral constraint value
    VectorXd lambda_normal(num_contacts); lambda_normal.setZero(); //  the normal constraint force
    VectorXd lambda_friction(2 * num_contacts); lambda_friction.setZero(); // the friction constraint force
    VectorXd lambda_bilateral(9 * num_bilateral_constraints); lambda_bilateral.setZero(); // the bilateral constraint force
    VectorXd lambda(lambda_bilateral.size() + lambda_normal.size() + lambda_friction.size()); lambda.setZero(); // the generalized constraint vector
    MatrixXd J_normal(num_contacts, num_elements * 3); J_normal.setZero(); // the normal jacobian matrix
    MatrixXd J_friction(2 * num_contacts, num_elements * 3); J_friction.setZero(); // the friction jacobian
    MatrixXd J_bilateral(9 * num_bilateral_constraints ,num_elements * 3); J_bilateral.setZero(); // the bilateral jacobian
    MatrixXd E(9*num_bilateral_constraints, 9*num_bilateral_constraints); E.setZero(); // the compliance matrix
    MatrixXd S(num_contacts, num_contacts); S.setZero(); // the normal compliance
    MatrixXd W(2 * num_contacts, 2 * num_contacts); W.setZero(); // the friction compliance

    VectorXd u(3 * num_elements); u.setZero() ; // the constrained velocity in generalized coordinates
    VectorXd u_tilde(3 * num_elements); u_tilde.setZero(); // the unconstrained velocity
    MatrixXd M_tilde = M;
    MatrixXd M_tilde_inv = M_inv;


    // fill the unconstrained positions
    int idxctr = 0;
    for (auto b: m_objects)
    {
        u_tilde.block(idxctr,0,3*b->m_numVerts,1) = b->xdot;
        idxctr += 3*b->m_numVerts;
    }

    ////////////////////////////////////////////////////////////////////////
    ////////////////// START THE NEWTON ITERATIONS ////////////////////////
    ////////////////////////////////////////////////////////////////////////

    for (int n = 0; n < m_settings.newtonIterations; n++)
    {
        cout << "Newton It = " << n << endl;

        for (auto b : m_objects)
        {
            int bidx = b->m_idx;
//            u.block<6,1>(6*bidx,0) = b->generalized_vel();
        }

        ////////////////////////////////////////////////////////////////////////////////////
        /////////////////////////// THE CONTACT CONSTRAINTS ////////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

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

        ////////////////////////////////////////////////////////////////////////////////////
        ///////////////////// THE SOFT KINEMATIC CONSTRAINTS //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

        for (auto obj : m_objects)
        {
            if (obj->type != DEFORMABLE) continue;
            auto& x = obj->x;
            for (auto c : obj->linear_hookean_constraints)
            {
                int tidx = c->tidx;
                VectorXd lambda_b = lambda_bilateral.block<6,1>(6*tidx,0);
                E.block<9,9>(tidx*9, tidx*9) = c->C;
                // the vertex indices that this constraint acts on
                int vidx[4]; vidx[0] = c->v0; vidx[1] = c->v1; vidx[2] = c->v2; vidx[3] = c->v3;
                // this constraint has 9 associated lagrange multipliers and acts on four vertices with 3 dof each
                // therefore the constraint jacobian is 9 x 12
                auto J_b = c->computeJacobian(x);
                for (int i = 0 ; i < 9; i++)
                {
                    for (int j = 0; j < 4; j++) {
                        J_bilateral.block<1, 3>(9 * tidx + i,3 * vidx[j]) = J_b.block<1,3>( i,3 * j);
                    }
                }

            }
        }

        cout << "Milestone " << endl;

        ////////////////////////////////////////////////////////////////////////////////////
        ///////////////////// THE HARD KINEMATIC CONSTRAINTS //////////////////////////////
        ///////////////////////////////////////////////////////////////////////////////////

        for (auto obj : m_objects)
        {
            if (obj->type != DEFORMABLE) continue;
            auto& x = obj->x;
            for (auto c : obj->fixed_constraints)
            {

            }
        }

        assert(!J_friction.hasNaN());
        assert(!J_normal.hasNaN());
        assert(!S.hasNaN());
        assert(!W.hasNaN());
        MatrixXd H = M_tilde;
        assert(!H.hasNaN());
        assert(J_normal.cols() == J_friction.cols());

        MatrixXd J(J_bilateral.rows() + J_normal.rows() + J_friction.rows(), J_normal.cols());
        J.setZero();
        J << J_bilateral , J_normal, J_friction;
        assert(!J.hasNaN());

        MatrixXd C(E.rows() + S.rows() + W.rows(),E.cols() + S.cols() + W.cols());
        C.setZero();
        C.block(0, 0, E.rows(), E.cols()) = E / pow(dt, 2);
        C.block(E.rows(), E.cols(), S.rows(), S.cols()) = S / pow(dt, 2);
        C.block( E.rows() + S.rows(), E.cols() + S.cols(), W.rows(), W.cols()) = W / dt;

        assert(!C.hasNaN());
        assert(!lambda.hasNaN());
        assert(J.rows() == lambda.size());
        assert(M_tilde.cols() == u.size());
        VectorXd g = M_tilde * (u - u_tilde) - dt * J.transpose() * lambda;
        assert(!g.hasNaN());
        MatrixXd I_b = MatrixXd::Identity(h_bilateral.size(),h_bilateral.size());
        I_b /= dt;
        MatrixXd I_n = MatrixXd::Identity(phi_normal.size(), phi_normal.size());
        I_n /= dt;
        MatrixXd I_f = MatrixXd::Identity(phi_friction.size(), phi_friction.size());
        MatrixXd I(I_b.rows() + I_n.rows() + I_f.rows(), I_b.cols() + I_n.cols() + I_f.cols());
        I.block(0, 0, I_b.rows(), I_b.cols()) = I_b;
        I.block(I_b.rows(), I_b.cols(), I_n.rows(), I_n.cols()) = I_n;
        I.block(I_b.rows() + I_n.rows(), I_b.cols()+ I_n.cols(), I_f.rows(), I_f.cols()) = I_f;
        VectorXd phi(h_bilateral.size() + phi_normal.size() + phi_friction.size());
        phi << h_bilateral , phi_normal, phi_friction;
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
        assert(!A.hasNaN());
        assert(J.cols() == H.rows() && H.cols() == g.size());
        assert(J.rows() == h.size());
        VectorXd b = (J * H_inv * g - h) / dt;
        assert(!b.hasNaN());

        // zero out any inactive contacts
        for (int cidx = 0; cidx < num_contacts; cidx++)
        {
            if (m_collisionDetector->m_contacts[cidx]->is_active() == 0 )
            {
                A.row(num_contacts + 2*cidx) = RowVectorXd::Zero(A.cols());
                A.row(num_contacts + 2*cidx + 1) = RowVectorXd::Zero(A.cols());
                b(num_contacts + 2*cidx) = 0;
                b(num_contacts + 2*cidx + 1) = 0;
            }
        }

        // solve this system with Gauss-Seidel
        VectorXd delta_lambda = A.fullPivLu().solve(b);

        cout << "DEBUG: Delta_lambda = " << delta_lambda.transpose() << endl;

        assert(!delta_lambda.hasNaN());
        auto delta_u = H_inv*(J.transpose()*delta_lambda*dt - g);
        assert(!delta_u.hasNaN());
        // perform a line search
        int t = 0.5;

        // update solution
        lambda += t*delta_lambda;
        u += t*delta_u;

        for (auto b : m_objects)
        {
            if (!b->is_static) {
                int bidx = b->m_idx;
                b->xdot = u;
                b->x += dt * b->xdot;
            }
            else
            {
                b->xdot.setZero();
            }
        }

        // set the constraint forces
        lambda_normal = lambda.block(0,0,lambda_normal.size(),1);
        lambda_friction = lambda.block(lambda_normal.size(),0,lambda_friction.size(),1);
    }
}