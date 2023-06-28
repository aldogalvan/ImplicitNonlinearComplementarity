#include "objects.hpp"

void RigidObject::set_local_pos(Eigen::Vector3d pos) {
    this->setLocalPos(pos);
    x_unconstrained = pos;
    x = pos;
}

MatrixXd RigidObject::kinematic_map_G()
{
    MatrixXd G(7,6); G.setZero();
    G.block<3,3>(0,0) = MatrixXd::Identity(3,3);
    G(3,3) =  -0.5*q.x();  G(3,4) =  -0.5*q.y(); G(3,5) =  -0.5*q.z();
    G(4,3) =  0.5*q.w();  G(4,4)  =   0.5*q.z(); G(4,5) =  -0.5*q.y();
    G(5,3) =  -0.5*q.z();  G(4,4) =  0.5*q.w(); G(4,5) =  0.5*q.x();
    G(6,3) =  0.5*q.y();  G(4,4)  =  -0.5*q.x(); G(4,5) =  0.5*q.w();
    return G;
}

MatrixXd RigidObject::kinematic_map_G(Quaterniond q)
{
    MatrixXd G(7,6); G.setZero();
    G.block<3,3>(0,0) = MatrixXd::Identity(3,3);
    G(3,3) =  -0.5*q.x();  G(3,4)  =  -0.5*q.y(); G(3,5) =  -0.5*q.z();
    G(4,3) =   0.5*q.w();  G(4,4)  =   0.5*q.z(); G(4,5) =  -0.5*q.y();
    G(5,3) =  -0.5*q.z();  G(4,4)  =   0.5*q.w(); G(4,5) =   0.5*q.x();
    G(6,3) =   0.5*q.y();  G(4,4)  =  -0.5*q.x(); G(4,5) =   0.5*q.w();
    return G;
}

MatrixXd RigidObject::kinematic_map_H()
{
    MatrixXd G(4,3); G.setZero();
    G(0,0) =  -0.5*q.x();   G(0,1) =  -0.5*q.y(); G(0,2) =  -0.5*q.z();
    G(1,0) =   0.5*q.w();   G(1,1) =   0.5*q.z(); G(1,2) =  -0.5*q.y();
    G(2,0) =  -0.5*q.z();   G(2,1) =   0.5*q.w(); G(2,2) =   0.5*q.x();
    G(3,0) =   0.5*q.y();   G(3,1) =  -0.5*q.x(); G(3,2) =   0.5*q.w();
    return G;
}

MatrixXd RigidObject::mass_matrix()
{
    MatrixXd M(6,6); M.setZero();
    if (!is_static)
    {
        M.block<3, 3>(0, 0) = mass * MatrixXd::Identity(3, 3);
        M.block<3, 3>(3, 3) = I;
    }

    return M;
}

MatrixXd RigidObject::inverse_mass_matrix()
{
    MatrixXd M_inv(6,6); M_inv.setZero();
    if (!is_static)
    {
        M_inv.block<3, 3>(0, 0) = (1 / mass) * MatrixXd::Identity(3, 3);
        M_inv.block<3, 3>(3, 3) = Iinv;
    }
    return M_inv;
}

void RigidObject::update_inertia_matrix()
{
    if( !is_static ){
        I = q * Ibody * q.inverse();
        Iinv = q * IbodyInv * q.inverse();}
    else{Iinv.setZero();}
}

void RigidObject::compute_inertia_matrix()
{
    Eigen::Matrix3d inertiaMatrix = Eigen::Matrix3d::Zero();

    for (int tidx = 0; tidx < triangles.rows(); tidx++)
    {
        Vector3i triangle = triangles.row(tidx);
        Eigen::Matrix3d triangleVertices;
        triangleVertices << vertices.row(triangle(0)),
                vertices.row(triangle(1)),
                vertices.row(triangle(2));

        Ibody = Eigen::Matrix3d::Zero();

        for (int i = 0; i < 3; i++)
        {
            Eigen::Vector3d vertex = triangleVertices.row(i);

            double x = vertex.x();
            double y = vertex.y();
            double z = vertex.z();

            double x2 = x * x;
            double y2 = y * y;
            double z2 = z * z;

            double xy = x * y;
            double xz = x * z;
            double yz = y * z;

            Ibody(0, 0) += y2 + z2;
            Ibody(1, 1) += x2 + z2;
            Ibody(2, 2) += x2 + y2;

            Ibody(0, 1) -= xy;
            Ibody(1, 0) -= xy;

            Ibody(0, 2) -= xz;
            Ibody(2, 0) -= xz;

            Ibody(1, 2) -= yz;
            Ibody(2, 1) -= yz;
        }

        inertiaMatrix += Ibody;
    }
    IbodyInv = Ibody.inverse();
}

void RigidObject::update_mesh_position()
{

    // Compute the transformed vertices
    MatrixXd transformedVertices = vertices;
    applyQuaternionRotation(transformedVertices,q_unconstrained);
    transformedVertices.rowwise() += x_unconstrained.transpose();

    // Update the vertices in the cMultiMesh object
    cMesh* mesh = this->getMesh(0);
    for (int j = 0; j < mesh->getNumVertices(); ++j) {
        mesh->m_vertices->setLocalPos(j, cVector3d(transformedVertices(j, 0),
                                                   transformedVertices(j, 1),
                                                   transformedVertices(j, 2)));
    }
}

void RigidObject::import_mesh_data()
{
    // Get the number of vertices and triangles in the mesh
    int numVertices = this->getNumVertices();
    int numTriangles = this->getNumTriangles();

    // Resize the Eigen matrices to accommodate the vertices and triangles
    vertices.resize(numVertices, 3);
    triangles.resize(numTriangles, 3);

    // Fill the vertices matrix
    for (int i = 0; i < numVertices; ++i) {

        const cVector3d& vertex = this->getVertexPos(i).eigen();
        vertices.row(i) << vertex.x(), vertex.y(), vertex.z();
    }

    // Fill the triangles matrix
    for (int i = 0; i < numTriangles; ++i) {
        cVector3d triangle (this->getMesh(0)->m_triangles->getVertexIndex0(i),
                            this->getMesh(0)->m_triangles->getVertexIndex1(i),
                            this->getMesh(0)->m_triangles->getVertexIndex2(i));
        triangles.row(i) << triangle.x(), triangle.y(), triangle.z();
    }

    // computes the inertia matrix
    compute_inertia_matrix();

}
