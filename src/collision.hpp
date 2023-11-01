
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP

#include <Eigen/Dense>
#include <vector>
#include "aabb.hpp"
#include "contact.hpp"
#include "objects.hpp"

using namespace Eigen;
using namespace std;


class CollisionDetector
{
public:

    CollisionDetector(vector<RigidObject*> objects)
    {
        m_objects = objects;
    }

    CollisionDetector(vector<DeformableObject*> objects)
    {

    }

    ~CollisionDetector();

    void computeCollisions(void);
    void clear(){for (auto c : m_contacts){delete c; m_contacts.clear(); }}
    void computeCollisionsRigidRigid(RigidObject* rigidObject1 , RigidObject* rigidObject2);
    void computeCollisionsRigidDeformable(RigidObject* rigidObject, DeformableObject* deformableObject);
    void computeCollisionsDeformableDeformable(DeformableObject* deformableObject1, DeformableObject* deformableObject2);
    static bool findCollisions(const Vector3d& object1_pos_start, const Vector3d& object1_pos_end, const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                        const Vector3d& object2_pos_start, const Vector3d& object2_pos_end, const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                        vector<Contact*>& collisions);

    static bool findCollisionsRigidRigid(const VectorXd& object1_start_position, const VectorXd& object1_end_position,
                               const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                               const VectorXd& object2_start_position, const VectorXd& object2_end_position,
                               const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                               vector<Contact*>& collisions, vector<cShapeLine*>& visualize_contact);
    static bool findCollisionsDeformableDeformable(const VectorXd& object1_start_vertices, const VectorXd& object1_end_vertices,
                                         const MatrixXi& object1_tris, const MatrixXi& object1_tets,
                                         const VectorXd& object2_start_vertices, const VectorXd& object2_end_vertices,
                                         const MatrixXi& object2_tris, const MatrixXi& object2_tets,
                                         vector<Contact*>& collisions);

    static bool findCollisionsRigidDeformable(VectorXd& object1_start_vertices, VectorXd& object1_end_vertices,
                                                   const MatrixXi& object1_tris, const MatrixXd& object2_start_vertices, const MatrixXd& object2_end_vertices,
                                                   const MatrixXi& object2_tris,
                                                   vector<Contact*>& collisions);

    static bool findCollisionsBroadPhase(const Eigen::VectorXd &object1_start_position,
                                         const Eigen::VectorXd &object1_end_position,
                                         const Eigen::MatrixXd &object1_vertices,
                                         const Eigen::MatrixXi &object1_tris,
                                         const Eigen::VectorXd &object2_start_position,
                                         const Eigen::VectorXd &object2_end_position,
                                         const Eigen::MatrixXd &object2_vertices,
                                         const Eigen::MatrixXi &object2_tris, std::vector<Collision>& potentialCollisions);

    // the objects in this detector
    vector<RigidObject*> m_objects;

    // the contacts detected by this detector
    vector<Contact*> m_contacts;
};



#endif //IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
