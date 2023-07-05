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
            m_objects.back()->set_local_pos(Vector3d(0,sqrt(3)/4*0.55,0.7));
            m_objects.back()->getMesh(0)->m_material->setBlue();
        }
        else if(i == 1)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj"));
            m_objects.back()->set_local_pos(Vector3d(-0.55/2,-sqrt(3)/4*0.55,0.7));
            m_objects.back()->getMesh(0)->m_material->setRed();
        }
        else if (i == 2)
        {
            m_objects.emplace_back(new RigidObject(i,"/home/agalvan-admin/ImplicitNonlinearComplementarity/resources/model.obj"));
            m_objects.back()->set_local_pos(Vector3d(0.55/2,-sqrt(3)/4*0.55,0.7));
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
        b->f = b->mass * Vector3d(0., 0, -.00981);
        b->tau.setZero();
        b->fc.setZero();
        b->tauc.setZero();
        b->m_contacts.clear();
    }

    // updates the inertia matrix
    for(auto b : m_objects)
    {
        b->R = b->q.toRotationMatrix();
        b->update_inertia_matrix();
    }

    // computes any collisions
    m_collisionDetector->computeCollisions();

    // computes the constraints

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

void RigidBodyDemo::computeConstraints()
{

}