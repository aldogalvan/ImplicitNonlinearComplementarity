#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_OBJECTS_HPP

#include <Eigen/Dense>
#include "chai3d.h"
#include "helper.hpp"

using namespace chai3d;
using namespace Eigen;

//struct GodObject
//{
//
//    GodObject()
//    {
//        q_eye.resize(6);
//        qdot_eye.resize(6);
//        qddot_eye.resize(6);
//        m.resize(6,6);
//        m.setIdentity();
//    }
//
//    void updateMeshPosition()
//    {
//        // Ensure the necessary data is not null
//        if (vis == nullptr || vertices == nullptr)
//            return;
//
//        // Compute the transformed vertices
//        Vector3d pos = q_eye.head(3);
//        Matrix3d rot = anglesToRotationMatrix(q_eye.tail(3));
//        MatrixXd transformedVertices = (rot) * (*vertices).transpose();
//        transformedVertices.colwise() += (pos);
//
//        // Update the vertices in the cMultiMesh object
//        cMesh* mesh = vis->getMesh(0);
//        for (int j = 0; j < mesh->getNumVertices(); ++j)
//        {
//            mesh->m_vertices->setLocalPos(j, cVector3d(transformedVertices(0, j),
//                                                       transformedVertices(1, j),
//                                                       transformedVertices(2, j)));
//        }
//    }
//
//    MatrixXd m;
//    VectorXd q_eye;
//    VectorXd q_eye_minus_one;
//    VectorXd qdot_eye;
//    VectorXd qddot_eye;
//
//    //
//    cMultiMesh* vis;
//    MatrixXi* triangles;
//    MatrixXd* vertices;
//    Vector3d* com;
//
//};

class RigidObject
{

public:

    RigidObject()
    {
        q_eye.resize(6);
        qdot_eye.resize(6);
        qddot_eye.resize(6);
        M.resize(3,3);
        M.setIdentity();
    }

    void updateMeshPosition()
    {
        // Ensure the necessary data is not null
        if (vis == nullptr || vertices == nullptr)
            return;

        // Compute the transformed vertices
        Vector3d pos = q_eye.head(3);
        Matrix3d rot = anglesToRotationMatrix(q_eye.tail(3));
        MatrixXd transformedVertices = (rot) * (*vertices).transpose();
        transformedVertices.colwise() += (pos);

        // Update the vertices in the cMultiMesh object
        cMesh* mesh = vis->getMesh(0);
        for (int j = 0; j < mesh->getNumVertices(); ++j) {
            mesh->m_vertices->setLocalPos(j, cVector3d(transformedVertices(0, j),
                                                       transformedVertices(1, j),
                                                       transformedVertices(2, j)));
        }
    }


    MatrixXd M;
    VectorXd q_eye;
    VectorXd q_eye_minus_one;
    VectorXd qdot_eye;
    VectorXd qddot_eye;
    //
    cMultiMesh* vis;
    MatrixXi* triangles;
    MatrixXd* vertices;
    Vector3d* com;

    bool is_static = false;
};

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_LCP_HPP
