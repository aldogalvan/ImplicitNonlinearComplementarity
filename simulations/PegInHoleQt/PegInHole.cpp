
#include "PegInHole.hpp"


void PegInHole::initialize()
{
    // create the meshes
    peg = new RigidObject();
    peg->newMesh();
    cCreateCylinder(peg->getMesh(0),1,1);
    peg->getMesh(0)->m_material->setBlue(); peg->getMesh(0)->m_material->setShininess(1);
    m_world->addChild(peg);

    // create the block
    block = new RigidObject();
    block->newMesh();
    cCreateBox(block->getMesh(0),1,1,1);
    block->getMesh(0)->m_material->setRed(); block->getMesh(0)->m_material->setShininess(1);
    m_world->addChild(block);
}

void PegInHole::step(double dt)
{

}

void PegInHole::solve_local()
{

}

void PegInHole::solve_global()
{

}

void PegInHole::update_haptics()
{

}

void PegInHole::update_graphics()
{

}