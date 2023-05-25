#include "collision.hpp"
#include "CTCD.h"

bool findCollisions(MatrixXd object1_start, MatrixXd object1_end, MatrixXi object1_tris,
                    MatrixXd object2_start, MatrixXd object2_end, MatrixXi object2_tris,
                    vector<ColInfo*>& collisions)
{
    bool flag = false;

    auto object1_aabb = buildAABB(object1_start, object1_end, object1_tris);
    auto object2_aabb = buildAABB(object2_start, object2_end, object2_tris);

    std::vector<Collision> potentialCollisions;
    intersect(object1_aabb,object2_aabb,potentialCollisions);

    for (const auto& it : potentialCollisions) {
        Eigen::Vector3i face = object2_tris.row(it.collidingTriangle2);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object1_tris(it.collidingTriangle1, vert);
            double t;

            Vector3d a0 = object1_start.row(vidx).transpose();
            Vector3d b0 = object2_start.row(face(0)).transpose();
            Vector3d c0 = object2_start.row(face(1)).transpose();
            Vector3d d0 = object2_start.row(face(2)).transpose();
            Vector3d a1 = object1_end.row(vidx).transpose();
            Vector3d b1 = object2_end.row(face(0)).transpose();
            Vector3d c1 = object2_end.row(face(1)).transpose();
            Vector3d d1 = object2_end.row(face(2)).transpose();

            // TODO Implement Vertex-Face CTCD
            if (CTCD::vertexFaceCTCD(
                    a0,b0,c0,d0,
                    a1,b1,c1,d1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new ColInfo);
                collisions.back()->t = t;
                Vector3d n0 = (c0 - b0).cross(d0 - b0);
                Vector3d n1 = (c1 - b1).cross(d1 - b1);
                Vector3d ncol = n0 + t*(n1 - n0);
                collisions.back()->normal = n1.normalized();
                collisions.back()->bodyIdxA = 0;
                collisions.back()->bodyIdxB = 1;
                flag = true;
            }
        }
    }

    for (const auto& it : potentialCollisions) {
        Eigen::Vector3i face = object1_tris.row(it.collidingTriangle1);
        for (int vert = 0; vert < 3; vert++) {
            int vidx = object2_tris(it.collidingTriangle2, vert);
            double t;

            Vector3d a0 = object2_start.row(vidx).transpose();
            Vector3d b0 = object1_start.row(face(0)).transpose();
            Vector3d c0 = object1_start.row(face(1)).transpose();
            Vector3d d0 = object1_start.row(face(2)).transpose();
            Vector3d a1 = object2_end.row(vidx).transpose();
            Vector3d b1 = object1_end.row(face(0)).transpose();
            Vector3d c1 = object1_end.row(face(1)).transpose();
            Vector3d d1 = object1_end.row(face(2)).transpose();

            if (CTCD::vertexFaceCTCD(
                    a0,b0,c0,d0,
                    a1,b1,c1,d1,
                    1e-6,
                    t))
            {
                collisions.emplace_back(new ColInfo);
                collisions.back()->t = t;
                Vector3d n0 = (c0 - b0).cross(d0 - b0);
                Vector3d n1 = (c1 - b1).cross(d1 - b1);
                Vector3d ncol = n0 + t*(n1 - n0);
                collisions.back()->normal = n1.normalized();
                collisions.back()->bodyIdxA = 1;
                collisions.back()->bodyIdxB = 0;
                flag = true;
            }
        }
    }
    return flag;
}