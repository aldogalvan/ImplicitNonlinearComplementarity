#include "objects.hpp"
#include "tetgen.h"


void RigidObject::set_local_pos(Eigen::Vector3d pos) {
    this->setLocalPos(pos);
    x_last = pos;
    x = pos;
}

void RigidObject::set_local_rot(Quaterniond rot)
{
    q_last = rot;
    q = rot;
}

MatrixXd RigidObject::kinematic_map_G()
{
    MatrixXd G(7,6); G.setZero();
    G.block<3,3>(0,0) = MatrixXd::Identity(3,3);
    G(3,3) =  -0.5*q.x();  G(3,4) =  -0.5*q.y(); G(3,5) =  -0.5*q.z();
    G(4,3) =  0.5*q.w();  G(4,4)  =   0.5*q.z(); G(4,5) =  -0.5*q.y();
    G(5,3) =  -0.5*q.z();  G(5,4) =  0.5*q.w(); G(5,5) =  0.5*q.x();
    G(6,3) =  0.5*q.y();  G(6,4)  =  -0.5*q.x(); G(6,5) =  0.5*q.w();
    return G;
}

MatrixXd RigidObject::kinematic_map_G(Quaterniond q)
{
    MatrixXd G(7,6); G.setZero();
    G.block<3,3>(0,0) = MatrixXd::Identity(3,3);
    G(3,3) =  -0.5*q.x();  G(3,4)  =  -0.5*q.y(); G(3,5) =  -0.5*q.z();
    G(4,3) =   0.5*q.w();  G(4,4)  =   0.5*q.z(); G(4,5) =  -0.5*q.y();
    G(5,3) =  -0.5*q.z();  G(5,4)  =   0.5*q.w(); G(5,5) =   0.5*q.x();
    G(6,3) =   0.5*q.y();  G(6,4)  =  -0.5*q.x(); G(6,5) =   0.5*q.w();
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
    Ibody = Eigen::Matrix3d::Zero();

    for (int tidx = 0; tidx < m_triangles.rows(); tidx++)
    {
        Vector3i triangle = m_triangles.row(tidx);
        Eigen::Matrix3d triangleVertices;
        triangleVertices << m_vertices.row(triangle(0)),
                m_vertices.row(triangle(1)),
                m_vertices.row(triangle(2));

        Eigen::Matrix3d inertiaMatrix = Eigen::Matrix3d::Zero();

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

            inertiaMatrix(0, 0) += y2 + z2;
            inertiaMatrix(1, 1) += x2 + z2;
            inertiaMatrix(2, 2) += x2 + y2;

            inertiaMatrix(0, 1) -= xy;
            inertiaMatrix(1, 0) -= xy;

            inertiaMatrix(0, 2) -= xz;
            inertiaMatrix(2, 0) -= xz;

            inertiaMatrix(1, 2) -= yz;
            inertiaMatrix(2, 1) -= yz;
        }

        Ibody += inertiaMatrix;
    }
    IbodyInv = Ibody.inverse();
}

void RigidObject::update_mesh_position()
{

    // Compute the transformed vertices
//    MatrixXd transformedVertices = m_vertices;
//    applyQuaternionRotation(transformedVertices,q);
//    transformedVertices.rowwise() += x.transpose();
//
//    // Update the vertices in the cMultiMesh object
//    cMesh* mesh = this->getMesh(0);
//    for (int j = 0; j < mesh->getNumVertices(); ++j) {
//        mesh->m_vertices->setLocalPos(j, cVector3d(transformedVertices(j, 0),
//                                                   transformedVertices(j, 1),
//                                                   transformedVertices(j, 2)));
//    }
    this->setLocalPos(x); this->setLocalRot(q.toRotationMatrix());
}

void RigidObject::import_mesh_data()
{
    // Get the number of vertices and triangles in the mesh
    int numVertices = this->getNumVertices();
    int numTriangles = this->getNumTriangles();

    // Resize the Eigen matrices to accommodate the vertices and triangles
    m_vertices.resize(numVertices, 3);
    m_triangles.resize(numTriangles, 3);

    // Fill the vertices matrix
    for (int i = 0; i < numVertices; ++i) {

        const cVector3d& vertex = this->getVertexPos(i).eigen();
        m_vertices.row(i) << vertex.x(), vertex.y(), vertex.z();
    }

    // Fill the triangles matrix
    for (int i = 0; i < numTriangles; ++i) {
        cVector3d triangle (this->getMesh(0)->m_triangles->getVertexIndex0(i),
                            this->getMesh(0)->m_triangles->getVertexIndex1(i),
                            this->getMesh(0)->m_triangles->getVertexIndex2(i));
        m_triangles.row(i) << triangle.x(), triangle.y(), triangle.z();
    }

    // computes the inertia matrix
    compute_inertia_matrix();

}


bool DeformableObject::create_tetrahedral_mesh(char* filename)
{

    tetgenio input;

    this->newMesh();
    cMesh* mesh = this->getMesh(0);

    // TetGen switches
    char TETGEN_SWITCHES[] = "pq1.414a0.002";

    if (input.load_off(filename))
    {
        // use TetGen to tetrahedralize our mesh
        tetgenio output;
        tetrahedralize(TETGEN_SWITCHES, &input, &output);

        m_vertices.resize(output.numberofpoints,3);
        m_numVerts = output.numberofpoints;

        // create a vertex in the object for each point of the result
        for (int p = 0, pi = 0; p < output.numberofpoints; ++p, pi += 3)
        {
            cVector3d point;
            point.x(output.pointlist[pi + 0]);
            point.y(output.pointlist[pi + 1]);
            point.z(output.pointlist[pi + 2]);

            m_vertices.row(p) = Eigen::RowVector3d(output.pointlist[pi + 0],
                                          output.pointlist[pi + 1],
                                          output.pointlist[pi + 2]);

            mesh->newVertex(point);
        }

        m_triangles.resize(output.numberoftrifaces,3);

        // create a triangle for each face on the surface
        for (int t = 0, ti = 0; t < output.numberoftrifaces; ++t, ti += 3)
        {
            cVector3d p[3];
            unsigned int vi[3];

            for (int i = 0; i < 3; ++i)
            {
                int tc = output.trifacelist[ti + i];
                vi[i] = tc;
                int pi = tc * 3;
                p[i].x(output.pointlist[pi + 0]);
                p[i].y(output.pointlist[pi + 1]);
                p[i].z(output.pointlist[pi + 2]);
            }

            m_triangles.row(t) = Eigen::RowVector3i(vi[0],vi[1],vi[2]);
            mesh->newTriangle(vi[0], vi[1], vi[2]);
        }

        m_tetrahedra.resize(output.numberoftetrahedra,4);
        m_numTets = output.numberoftetrahedra;

        for (int t = 0, ti = 0; t < output.numberoftetrahedra; ++t, ti += 4)
        {
            int v0 = output.tetrahedronlist[ti + 0];
            int v1 = output.tetrahedronlist[ti + 1];
            int v2 = output.tetrahedronlist[ti + 2];
            int v3 = output.tetrahedronlist[ti + 3];

            Eigen::RowVector4i tetrahedron;
            tetrahedron[0] = v0;
            tetrahedron[1] = v1;
            tetrahedron[2] = v2;
            tetrahedron[3] = v3;

            m_tetrahedra.row(t) = (tetrahedron);
        }
        return true;
    }
    else {
        cout << "FAILED TO LOAD TETRAHEDRAL MESH" << endl;
        return false;
    }

    // compute the rest volume
    compute_rest_volumes();
    compute_rest_deformation_gradient();
    compute_elasticity_matrix();
}

void DeformableObject::compute_rest_volumes()
{
    int numTets = m_tetrahedra.rows();
    volume_0.resize(numTets);
    for (int tidx = 0; tidx < numTets; tidx++)
    {
        int t0= m_tetrahedra(tidx,0); int t1= m_tetrahedra(tidx,1); int t2= m_tetrahedra(tidx,2); int t3= m_tetrahedra(tidx,3);
        Vector3d p0 = m_vertices.row(t0); Vector3d p1 = m_vertices.row(t1); Vector3d p2 = m_vertices.row(t2); Vector3d p3 = m_vertices.row(t3);
        volume_0(tidx) = computeTetrahedronVolume(p0,p1,p2,p3);
    }
}

void DeformableObject::compute_rest_deformation_gradient()
{
    int numTets = m_tetrahedra.rows();
    gradient_0.resize(3*numTets,3);
    for(int tidx = 0; tidx < numTets; tidx++)
    {
        int t0 = m_tetrahedra(tidx,0); int t1 = m_tetrahedra(tidx,1); int t2 = m_tetrahedra(tidx,2); int t3 = m_tetrahedra(tidx,3);
        Vector3d v0 = m_vertices.row(t0);
        Vector3d v1 = m_vertices.row(t1);
        Vector3d v2 = m_vertices.row(t2);
        Vector3d v3 = m_vertices.row(t3);
        MatrixXd r_(3,3);
        r_.col(0) = v0 - v1;
        r_.col(1) = v0 - v2;
        r_.col(2) = v0 - v3;
        gradient_0.block<3,3>(3*tidx,0) = r_;
    }
}

void DeformableObject::set_material_stiffness_matrix(double E, double nu)
{
    K.resize(6,6);

    // Calculate material constants
    double lambda = E * nu / ((1 + nu) * (1 - 2 * nu));
    double mu = E / (2 * (1 + nu));

    // Create the material stiffness matrix
    MatrixXd materialMatrix(6, 6);
    materialMatrix << lambda + 2 * mu, lambda, lambda, 0, 0, 0,
            lambda, lambda + 2 * mu, lambda, 0, 0, 0,
            lambda, lambda, lambda + 2 * mu, 0, 0, 0,
            0, 0, 0, mu, 0, 0,
            0, 0, 0, 0, mu, 0,
            0, 0, 0, 0, 0, mu;

}

void DeformableObject::compute_elasticity_matrix()
{
    int numTet = m_tetrahedra.rows();
    E.resize(numTet*6,numTet*6);
    auto K_inv = K.inverse();
    for (int tidx = 0; tidx < numTet; tidx++)
    {
        int v0 = m_tetrahedra(tidx,0);
        int v1 = m_tetrahedra(tidx,1);
        int v2 = m_tetrahedra(tidx,2);
        int v3 = m_tetrahedra(tidx,3);

        Vector3d p0 = m_vertices.row(v0);
        Vector3d p1 = m_vertices.row(v1);
        Vector3d p2 = m_vertices.row(v2);
        Vector3d p3 = m_vertices.row(v3);

        double vol = volume_0(tidx);
        E.block<6,6>(6*tidx,6*tidx) = vol*K_inv;
    }
}

void DeformableObject::update_mesh_position()
{
    cMesh* mesh = this->getMesh(0);
    int numVerts = m_vertices.rows();
    for (int vidx = 0; vidx < numVerts; vidx++)
    {
        mesh->m_vertices->setLocalPos(vidx,m_vertices(vidx,0),m_vertices(vidx,1),m_vertices(vidx,2));
    }
}

MatrixXd DeformableObject::mass_matrix()
{
    return mass*MatrixXd::Identity(3*m_numVerts,3*m_numVerts) / m_numVerts;
}

MatrixXd DeformableObject::inverse_mass_matrix()
{
    return m_numVerts*MatrixXd::Identity(3*m_numVerts,3*m_numVerts) / mass;
}