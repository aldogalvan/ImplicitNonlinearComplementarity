#include "collision.hpp"
#include "CTCD.h"
#include <iostream>
#include "helper.hpp"



inline Vector3d refinePoint(const Vector3d& p, const Vector3d& objectPosition, const Quaterniond& objectRotation)
{
    // Create the object's transformation matrix
    Affine3d objectTransform;
    objectTransform.translation() = objectPosition;
    objectTransform.linear() = objectRotation.toRotationMatrix();

    // Apply the inverse transformation to the point
    Vector3d refinedPoint = objectTransform.inverse() * p;

    return refinedPoint;
}

Vector3d closestPointOnTriangle(const Vector3d& A, const Vector3d& B, const Vector3d& C, const Vector3d& p) {
    Vector3d AB = B - A;
    Vector3d AC = C - A;
    Vector3d AP = p - A;

    double dABAB = AB.dot(AB);
    double dABAC = AB.dot(AC);
    double dACAC = AC.dot(AC);
    double dAPAB = AP.dot(AB);
    double dAPAC = AP.dot(AC);

    double denom = dABAB * dACAC - dABAC * dABAC;

    double s = (dABAB * dAPAC - dABAC * dAPAB) / denom;
    double t = (dACAC * dAPAB - dABAC * dAPAC) / denom;

    if (s < 0.0) {
        s = 0.0;
        t = dAPAC / dACAC;
    } else if (t < 0.0) {
        t = 0.0;
        s = dAPAB / dABAB;
    } else if (s + t > 1.0) {
        s = 1.0 - t;
    }

    return A + s * AB + t * AC;
}

void CollisionDetector::computeCollisions()
{
    // delete the current contacts
    for (auto c : m_contacts) {delete c;}
    m_contacts.clear();

    // Nested loops to iterate through the vector
    for (size_t i = 0; i < m_objects.size(); ++i) {
        for (size_t j = i + 1; j < m_objects.size(); ++j)
        {
            assert(m_objects[i] != NULL && m_objects[j] != NULL);
            // Call the function for each pair of objects
            if (m_objects[i]->type == RIGID && m_objects[j]->type == RIGID)
                computeCollisionsRigidRigid(dynamic_cast<RigidObject*>(m_objects[i]),dynamic_cast<RigidObject*>(m_objects[j]));
            else if (m_objects[i]->type == DEFORMABLE && m_objects[j]->type == DEFORMABLE)
                computeCollisionsDeformableDeformable(dynamic_cast<DeformableObject*>(m_objects[i]),dynamic_cast<DeformableObject*>(m_objects[j]));
            else if (m_objects[i]->type == RIGID && m_objects[j]->type == DEFORMABLE)
                computeCollisionsRigidDeformable(dynamic_cast<RigidObject*>(m_objects[i]),dynamic_cast<DeformableObject*>(m_objects[j]));
            else if (m_objects[i]->type == DEFORMABLE && m_objects[j]->type == RIGID)
                computeCollisionsRigidDeformable(dynamic_cast<RigidObject*>(m_objects[j]),dynamic_cast<DeformableObject*>(m_objects[i]));
        }
    }
}

void CollisionDetector::computeCollisionsRigidRigid(RigidObject* object1, RigidObject* object2)
{
    // declare Object 1
    const MatrixXd& object1_vertices = object1->vertices();
    const MatrixXi& object1_triangles = object1->triangles();
    Vector3d& object1_pos_start = object1->x_last; Vector3d& object1_pos_end = object1->x;
    Quaterniond& object1_rot_start = object1->q_last; Quaterniond& object1_rot_end = object1->q;

    // declare Object 2
    const MatrixXd& object2_vertices = object2->vertices();
    const MatrixXi& object2_triangles = object2->triangles();
    Vector3d& object2_pos_start = object2->x_last; Vector3d& object2_pos_end = object2->x;
    Quaterniond& object2_rot_start = object2->q_last; Quaterniond& object2_rot_end = object2->q;

    // object1 vertices
    MatrixXd object1_start, object1_end;
    object1_start = object1_vertices; object1_end = object1_vertices;
    applyQuaternionRotation(object1_start,object1_rot_start);
    applyQuaternionRotation(object1_end,object1_rot_end);
    object1_start.rowwise() += object1_pos_start.transpose();
    object1_end.rowwise() += object1_pos_end.transpose();

    // object2 vertices
    MatrixXd object2_start, object2_end;
    object2_start = object2_vertices; object2_end = object2_vertices;
    applyQuaternionRotation(object2_start,object2_rot_start);
    applyQuaternionRotation(object2_end,object2_rot_end);
    object2_start.rowwise() += object2_pos_start.transpose();
    object2_end.rowwise() += object2_pos_end.transpose();

    auto object1_aabb = buildAABB(object1_start, object1_end, object1_triangles);
    auto object2_aabb = buildAABB(object2_start, object2_end, object2_triangles);

    std::vector<Collision> potentialCollisions;
    intersect(object1_aabb,object2_aabb,potentialCollisions);


    //////////////////////////////////////////////////////////
    ///// HERE OBJECT 2 OWNS THE TRIANGLE ////////////////////
    //////////////////////////////////////////////////////////

    for (const auto& it : potentialCollisions) {
        Vector3i face = object2_triangles.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object1_triangles(it.collidingTriangle1, vert);
            double t;

            Vector3d v0 = object1_start.row(vidx).transpose();
            Vector3d a0 = object2_start.row(face(0)).transpose();
            Vector3d b0 = object2_start.row(face(1)).transpose();
            Vector3d c0 = object2_start.row(face(2)).transpose();
            Vector3d v1 = object1_end.row(vidx).transpose();
            Vector3d a1 = object2_end.row(face(0)).transpose();
            Vector3d b1 = object2_end.row(face(1)).transpose();
            Vector3d c1 = object2_end.row(face(2)).transpose();

            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                m_contacts.emplace_back(new Contact(object1,object2));
                m_contacts.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                m_contacts.back()->n = n1;
                m_contacts.back()->objectIdxA = object1->m_idx;
                m_contacts.back()->objectIdxB = object2->m_idx;
                m_contacts.back()->depth = (v1 - a1).dot(n1);

                // find the contacting point
                Vector3d P = v1 - m_contacts.back()->depth * n1;
                m_contacts.back()->contact_pt = P;
                m_contacts.back()->contact_wrt_objectA = refinePoint(v1,object1_pos_end,object1_rot_end);
                m_contacts.back()->contact_wrt_objectB = refinePoint(P,object2_pos_end,object2_rot_end);
                m_contacts.back()->compute_contact_basis();

            }
        }
    }

    //////////////////////////////////////////////////////////
    ///// HERE OBJECT 1 OWNS THE TRIANGLE ////////////////////
    //////////////////////////////////////////////////////////

    for (const auto& it : potentialCollisions) {
        Vector3i face = object1_triangles.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object2_triangles(it.collidingTriangle2, vert);
            double t;
            Vector3d v0 = object2_start.row(vidx).transpose();
            Vector3d a0 = object1_start.row(face(0)).transpose();
            Vector3d b0 = object1_start.row(face(1)).transpose();
            Vector3d c0 = object1_start.row(face(2)).transpose();
            Vector3d v1 = object2_end.row(vidx).transpose();
            Vector3d a1 = object1_end.row(face(0)).transpose();
            Vector3d b1 = object1_end.row(face(1)).transpose();
            Vector3d c1 = object1_end.row(face(2)).transpose();

            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {

                m_contacts.emplace_back(new Contact(object2,object1));
                m_contacts.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                m_contacts.back()->n = n1;
                m_contacts.back()->objectIdxA = object2->m_idx;
                m_contacts.back()->objectIdxB = object1->m_idx;
                m_contacts.back()->depth = ((v1 - a1).dot(n1));
                // find the contacting point
                Vector3d P = v1 - m_contacts.back()->depth * n1;
                m_contacts.back()->contact_pt = P;
                m_contacts.back()->contact_wrt_objectA = refinePoint(v1,object2_pos_end,object2_rot_end);
                m_contacts.back()->contact_wrt_objectB = refinePoint(P,object1_pos_end,object1_rot_end);
                m_contacts.back()->compute_contact_basis();

            }
        }
    }
}

void CollisionDetector::computeCollisionsRigidDeformable(RigidObject* rigidObject, DeformableObject* deformableObject){cout << "NOT WORKING" << endl;}
void CollisionDetector::computeCollisionsDeformableDeformable(DeformableObject* deformableObject1, DeformableObject* deformableObject2){cout << "NOT WORKING" << endl;}

bool CollisionDetector::findCollisions(const Vector3d& object1_pos_start, const Vector3d& object1_pos_end, const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                    const Vector3d& object2_pos_start, const Vector3d& object2_pos_end, const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                    vector<Contact*>& collisions)
{
    bool flag = false;

    // object1 vertices
    MatrixXd object1_start, object1_end;
    object1_start = object1_vertices; object1_end = object1_vertices;
    object1_start.rowwise() += object1_pos_start.transpose();
    object1_end.rowwise() += object1_pos_end.transpose();

    // object2 vertices
    MatrixXd object2_start, object2_end;
    object2_start = object2_vertices; object2_end = object2_vertices;
    object2_start.rowwise() += object2_pos_start.transpose();
    object2_end.rowwise() += object2_pos_end.transpose();

    auto object1_aabb = buildAABB(object1_start, object1_end, object1_tris);
    auto object2_aabb = buildAABB(object2_start, object2_end, object2_tris);

    std::vector<Collision> potentialCollisions;
    intersect(object1_aabb,object2_aabb,potentialCollisions);

    for (const auto& it : potentialCollisions) {
        Vector3i face = object2_tris.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object1_tris(it.collidingTriangle1, vert);
            double t;


            Vector3d v0 = object1_start.row(vidx).transpose();
            Vector3d a0 = object2_start.row(face(0)).transpose();
            Vector3d b0 = object2_start.row(face(1)).transpose();
            Vector3d c0 = object2_start.row(face(2)).transpose();
            Vector3d v1 = object1_end.row(vidx).transpose();
            Vector3d a1 = object2_end.row(face(0)).transpose();
            Vector3d b1 = object2_end.row(face(1)).transpose();
            Vector3d c1 = object2_end.row(face(2)).transpose();

            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = 0;
                collisions.back()->objectIdxB = -1;
                collisions.back()->depth = (v1 - a1).dot(n1);
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                collisions.back()->contact_wrt_objectA = v1 - object1_pos_end;
                collisions.back()->contact_wrt_objectB = P - object2_pos_end;
                collisions.back()->compute_contact_basis();
                flag = true;
            }
        }
    }

    for (const auto& it : potentialCollisions) {
        Vector3i face = object1_tris.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object2_tris(it.collidingTriangle2, vert);
            double t;
            Vector3d v0 = object2_start.row(vidx).transpose();
            Vector3d a0 = object1_start.row(face(0)).transpose();
            Vector3d b0 = object1_start.row(face(1)).transpose();
            Vector3d c0 = object1_start.row(face(2)).transpose();
            Vector3d v1 = object2_end.row(vidx).transpose();
            Vector3d a1 = object1_end.row(face(0)).transpose();
            Vector3d b1 = object1_end.row(face(1)).transpose();
            Vector3d c1 = object1_end.row(face(2)).transpose();

            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = -1;
                collisions.back()->objectIdxB = 0;
                collisions.back()->depth = ((v1 - a1).dot(n1));
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                collisions.back()->contact_wrt_objectA = v1 - object2_pos_end;
                collisions.back()->contact_wrt_objectB = P - object1_pos_end;
                collisions.back()->compute_contact_basis();
                flag = true;
            }
        }
    }
    return flag;
}

bool CollisionDetector::findCollisionsRigidRigid(const VectorXd& object1_start_position, const VectorXd& object1_end_position,
                    const MatrixXd& object1_vertices, const MatrixXi& object1_tris,
                   const VectorXd& object2_start_position, const VectorXd& object2_end_position,
                    const MatrixXd& object2_vertices, const MatrixXi& object2_tris,
                    vector<Contact*>& collisions)
{
    bool flag = false;

    // object1 vertices
    VectorXd object1_pos_start = object1_start_position.block<3,1>(0,0);
    VectorXd object1_pos_end = object1_end_position.block<3,1>(0,0);
    Quaterniond object1_rot_start(object1_start_position(3),object1_start_position(4),object1_start_position(5),object1_start_position(6));
    Quaterniond object1_rot_end(object1_end_position(3),object1_end_position(4),object1_end_position(5),object1_end_position(6));
    MatrixXd object1_start, object1_end;
    object1_start = object1_vertices; object1_end = object1_vertices;
    applyQuaternionRotation(object1_start,object1_rot_start);
    applyQuaternionRotation(object1_end,object1_rot_end);
    object1_start.rowwise() += object1_pos_start.transpose();
    object1_end.rowwise() += object1_pos_end.transpose();

    // object2 vertices
    VectorXd object2_pos_start = object2_start_position.block<3,1>(0,0);
    VectorXd object2_pos_end = object2_end_position.block<3,1>(0,0);
    Quaterniond object2_rot_start(object2_start_position(3),object2_start_position(4),object2_start_position(5),object2_start_position(6));
    Quaterniond object2_rot_end(object2_end_position(3),object2_end_position(4),object2_end_position(5),object2_end_position(6));
    MatrixXd object2_start, object2_end;
    object2_start = object2_vertices; object2_end = object2_vertices;
    applyQuaternionRotation(object2_start,object2_rot_start);
    applyQuaternionRotation(object2_end,object2_rot_end);
    object2_start.rowwise() += object2_pos_start.transpose();
    object2_end.rowwise() += object2_pos_end.transpose();

    auto object1_aabb = buildAABB(object1_start, object1_end, object1_tris);
    auto object2_aabb = buildAABB(object2_start, object2_end, object2_tris);

    std::vector<Collision> potentialCollisions;
    intersect(object1_aabb,object2_aabb,potentialCollisions);

    for (const auto& it : potentialCollisions) {
        Vector3i face = object2_tris.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object1_tris(it.collidingTriangle1, vert);
            double t;

            Vector3d v0 = object1_start.row(vidx).transpose();
            Vector3d a0 = object2_start.row(face(0)).transpose();
            Vector3d b0 = object2_start.row(face(1)).transpose();
            Vector3d c0 = object2_start.row(face(2)).transpose();
            Vector3d v1 = object1_end.row(vidx).transpose();
            Vector3d a1 = object2_end.row(face(0)).transpose();
            Vector3d b1 = object2_end.row(face(1)).transpose();
            Vector3d c1 = object2_end.row(face(2)).transpose();

            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = 0;
                collisions.back()->objectIdxB = -1;
                collisions.back()->depth = (v1 - a1).dot(n1);
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                collisions.back()->contact_wrt_objectA = refinePoint(v1,object1_pos_end,object1_rot_end);//v1 - object1_pos_end;
                collisions.back()->contact_wrt_objectB = refinePoint(P,object2_pos_end,object2_rot_end);//P - object2_pos_end;
                flag = true;
            }
        }
    }

    for (const auto& it : potentialCollisions) {
        Vector3i face = object1_tris.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object2_tris(it.collidingTriangle2, vert);
            double t;
            Vector3d v0 = object2_start.row(vidx).transpose();
            Vector3d a0 = object1_start.row(face(0)).transpose();
            Vector3d b0 = object1_start.row(face(1)).transpose();
            Vector3d c0 = object1_start.row(face(2)).transpose();
            Vector3d v1 = object2_end.row(vidx).transpose();
            Vector3d a1 = object1_end.row(face(0)).transpose();
            Vector3d b1 = object1_end.row(face(1)).transpose();
            Vector3d c1 = object1_end.row(face(2)).transpose();

            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = -1;
                collisions.back()->objectIdxB = 0;
                collisions.back()->depth = ((v1 - a1).dot(n1));
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                collisions.back()->contact_wrt_objectA = refinePoint(v1,object2_pos_end,object2_rot_end);//v1 - object2_pos_end;
                collisions.back()->contact_wrt_objectB = refinePoint(P,object1_pos_end,object1_rot_end);//P - object1_pos_end;
                flag = true;
            }
        }
    }
    return flag;
}

bool CollisionDetector::findCollisionsDeformableDeformable(const VectorXd& object1_start_vertices, const VectorXd& object1_end_vertices,
                                             const MatrixXi& object1_tris, const MatrixXi& object1_tets,
                                            const VectorXd& object2_start_vertices, const VectorXd& object2_end_vertices,
                                            const MatrixXi& object2_tris, const MatrixXi& object2_tets,
                                            vector<Contact*>& collisions)
{
    bool flag = false;

    // object1 vertices

    const MatrixXd& object1_start = object1_start_vertices; const MatrixXd& object1_end = object1_end_vertices;
    // object2 vertices
    const MatrixXd& object2_start = object2_start_vertices; const MatrixXd& object2_end = object2_end_vertices;

    auto object1_aabb = buildAABB(object1_start, object1_end, object1_tris);
    auto object2_aabb = buildAABB(object2_start, object2_end, object2_tris);

    std::vector<Collision> potentialCollisions;
    intersect(object1_aabb,object2_aabb,potentialCollisions);

    for (const auto& it : potentialCollisions) {
        Vector3i face = object2_tris.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object1_tris(it.collidingTriangle1, vert);
            double t;

            Vector3d v0 = object1_start.row(vidx).transpose();
            Vector3d a0 = object2_start.row(face(0)).transpose();
            Vector3d b0 = object2_start.row(face(1)).transpose();
            Vector3d c0 = object2_start.row(face(2)).transpose();
            Vector3d v1 = object1_end.row(vidx).transpose();
            Vector3d a1 = object2_end.row(face(0)).transpose();
            Vector3d b1 = object2_end.row(face(1)).transpose();
            Vector3d c1 = object2_end.row(face(2)).transpose();

            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = 0;
                collisions.back()->objectIdxB = -1;
                collisions.back()->depth = (v1 - a1).dot(n1);
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                flag = true;
            }
        }
    }

    for (const auto& it : potentialCollisions) {
        Vector3i face = object1_tris.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object2_tris(it.collidingTriangle2, vert);
            double t;
            Vector3d v0 = object2_start.row(vidx).transpose();
            Vector3d a0 = object1_start.row(face(0)).transpose();
            Vector3d b0 = object1_start.row(face(1)).transpose();
            Vector3d c0 = object1_start.row(face(2)).transpose();
            Vector3d v1 = object2_end.row(vidx).transpose();
            Vector3d a1 = object1_end.row(face(0)).transpose();
            Vector3d b1 = object1_end.row(face(1)).transpose();
            Vector3d c1 = object1_end.row(face(2)).transpose();

            if (CTCD::vertexFaceCTCD(
                    v0,a0,b0,c0,
                    v1,a1,b1,c1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new Contact);
                collisions.back()->t = t;
                Vector3d n1 = ((b1 - a1).cross(c1 - a1)).normalized();
                collisions.back()->n = n1;
                collisions.back()->objectIdxA = -1;
                collisions.back()->objectIdxB = 0;
                collisions.back()->depth = ((v1 - a1).dot(n1));
                // find the contacting point
                Vector3d P = v1 - collisions.back()->depth * n1;
                collisions.back()->contact_pt = P;
                flag = true;
            }
        }
    }
    return flag;
}