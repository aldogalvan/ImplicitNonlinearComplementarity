
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP

#include <Eigen/Dense>
#include <vector>
#include "objects.hpp"
#include "aabb.hpp"
#include "contact.hpp"

using namespace Eigen;
using namespace std;

class CollisionDetector
{
public:

    CollisionDetector(vector<RigidObject*> objects)
    {
        m_objects = objects;
    }

    ~CollisionDetector();

    void computeCollisions(void);
    void computeCollisions(RigidObject* object1 , RigidObject* object2);
    static bool findCollisions(const Vector3d& object1_pos_start, const Vector3d& object1_pos_end, const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                        const Vector3d& object2_pos_start, const Vector3d& object2_pos_end, const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                        vector<Contact*>& collisions);

    static bool findCollisionsRigid(const VectorXd& object1_start_position, const VectorXd& object1_end_position,
                               const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                               const VectorXd& object2_start_position, const VectorXd& object2_end_position,
                               const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                               vector<Contact*>& collisions);

    static bool findCollisionsDeformable(const VectorXd& object1_start_vertices, const VectorXd& object1_end_vertices,
                                         const MatrixXi& object1_tris, const MatrixXi& object1_tets,
                                         const VectorXd& object2_start_vertices, const VectorXd& object2_end_vertices,
                                         const MatrixXi& object2_tris, const MatrixXi& object2_tets,
                                         vector<Contact*>& collisions);

    // the objects in this detector
    vector<RigidObject*> m_objects;

    // the contacts detected by this detector
    vector<Contact*> m_contacts;
};



#endif //IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
