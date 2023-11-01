#include "objects.hpp"
#include "tetgen.h"
#include <igl/readOBJ.h>
#include <igl/readSTL.h>
#include <igl/readOFF.h>
#include <igl/centroid.h>
#include <random>


void generateRandomPastelColor(double& red, double& green, double& blue) {
    // Create a random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> distribution(0.5, 1.0);  // Limit to light shades in the [0.7, 1.0] range

    // Generate random values for the RGB components
    red = distribution(gen);
    green = distribution(gen);
    blue = distribution(gen);
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// BASE OBJECT /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void Object::loadMesh(const std::string& filename)
{
    // Determine the file extension
    size_t dotPos = filename.find_last_of('.');
    if (dotPos == std::string::npos)
    {
        std::cerr << "Invalid file: " << filename << std::endl;
        return;
    }

    std::string extension = filename.substr(dotPos + 1);
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    // Load the mesh based on the file extension
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd N;

    if (extension == "obj")
    {
        if (!igl::readOBJ(filename, V, F))
        {
            std::cerr << "Error loading OBJ file: " << filename << std::endl;
            return;
        }
    }
//    else if (extension == "stl")
//    {
//        if (!igl::readSTL(filename, V, F, N))
//        {
//            std::cerr << "Error loading STL file: " << filename << std::endl;
//            return;
//        }
//    }
    else if (extension == "off")
    {
        if (!igl::readOFF(filename, V, F))
        {
            std::cerr << "Error loading OFF file: " << filename << std::endl;
            return;
        }
    }
    else
    {
        std::cerr << "Invalid file type: " << extension << std::endl;
        return;
    }

    // Calculate the centroid (mean) of the vertices
    Eigen::Vector3d centroid;
    for (int k = 0 ; k < V.rows(); k++)
    {
        centroid += V.row(k);
    }
    centroid /= V.rows();

    // Translate the vertices to center the mesh around the origin
    V.rowwise() -= centroid.transpose();

    // Store the loaded vertex and triangle information in member variables
    m_vertices = V;
    m_triangles = F;
    m_normals = N;

    // Optionally, compute normals if needed
    // igl::per_vertex_normals(V, F, m_normals);

    // transfers the data to chai3d
    initializeVisualizer();
}

void Object::scaleMesh(const double s)
{
    m_vertices *= s;
    this->scale(s);
}


void Object::initializeVisualizer()
{
    // create a new mesh
    this->newMesh();
    cMesh* mesh = getMesh(0);

    // create triangles and vertices for the mesh
    for (int vidx = 0; vidx < m_vertices.rows(); vidx++)
    {
        mesh->newVertex(cVector3d(m_vertices.row(vidx)));
    }

    for (int tidx = 0; tidx < m_triangles.rows(); tidx++)
    {
        mesh->newTriangle(m_triangles(tidx,0),m_triangles(tidx,1),m_triangles(tidx,2));
    }
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// RIGID OBJECT /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

VectorXd RigidObject::get_configuration() {
    VectorXd ret(7);
    ret.head(3) = x;
    ret(3) = q.w();
    ret(4) = q.x();
    ret(5) = q.y();
    ret(6) = q.z();
    return ret;
}

VectorXd RigidObject::get_configuration_unconstrained() {
    VectorXd ret(7);
    ret.head(3) = x_tilde;
    ret(3) = q_tilde.w();
    ret(4) = q_tilde.x();
    ret(5) = q_tilde.y();
    ret(6) = q_tilde.z();    return ret;
}

void RigidObject::set_local_pos(Eigen::Vector3d pos) {
    this->setLocalPos(pos);
    x_tilde = pos;
    x = pos;
}

void RigidObject::set_local_rot(Quaterniond rot)
{
    q_tilde = rot;
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

VectorXd RigidObject::generalized_pos()
{
    q.normalize();
    VectorXd ret;
    ret.head(3) = x;
}

VectorXd RigidObject::generalized_vel()
{
    VectorXd u(6);
    u.head(3) = xdot;
    u.tail(3) = omega;
    return u;
}

MatrixXd RigidObject::generalized_mass()
{
    MatrixXd M(6,6); M.setZero();
    M.block<3, 3>(0, 0) = mass * MatrixXd::Identity(3, 3);
    M.block<3, 3>(3, 3) = I;


    return M;
}

MatrixXd RigidObject::generalized_mass_inverse()
{
    MatrixXd M_inv(6,6); M_inv.setZero();
    M_inv.block<3, 3>(0, 0) = (1 / mass) * MatrixXd::Identity(3, 3);
    M_inv.block<3, 3>(3, 3) = Iinv;

    return M_inv;
}

void RigidObject::update_inertia_matrix()
{
    I = q * Ibody * q.inverse();
    Iinv = q * IbodyInv * q.inverse();
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
    this->setLocalPos(x); this->setLocalRot(q.toRotationMatrix());
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// DEFORMABLE OBJECT /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

bool DeformableObject::create_tetrahedral_mesh(char* filename)
{
    tetgenio input;

    this->newMesh();
    cMesh* mesh = this->getMesh(0);

    // TetGen switches
    char TETGEN_SWITCHES[] = "pq1.414a0.2";

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

        // set the vertex positions
        x = flattenMatrix(m_vertices); x_tilde = x;
        xdot = VectorXd::Zero(m_numVerts*3);
        xdot_tilde = VectorXd::Zero(m_numVerts*3);
        f = VectorXd::Zero(m_numVerts*3);
        fc = VectorXd::Zero(m_numVerts*3);

        // creates the mass matrices

        M = MatrixXd::Zero(m_numVerts*3,m_numVerts*3);
        M_inv = MatrixXd::Zero(m_numVerts*3,m_numVerts*3);
        M_sparse = SparseMatrix<double>(m_numVerts*3,m_numVerts*3);
        M_inv_sparse = SparseMatrix<double>(m_numVerts*3,m_numVerts*3);
        for (int i = 0 ; i < 3*m_numVerts; i++)
        {
            M_sparse.insert(i,i) = mass / m_numVerts;
            M_inv_sparse.insert(i,i) = m_numVerts / mass;
            M(i,i) = mass/m_numVerts;
            M_inv(i,i) = m_numVerts / mass;
        }


        return true;
    }
    else {
        cout << "FAILED TO LOAD TETRAHEDRAL MESH" << endl;
        return false;
    }
}

void DeformableObject::update_mesh_position()
{
    cMesh* mesh = this->getMesh(0);
    int numVerts = m_vertices.rows();
    for (int vidx = 0; vidx < numVerts; vidx++)
    {
        mesh->m_vertices->setLocalPos(vidx,x(3*vidx),x(3*vidx + 1),x(3*vidx + 2));
    }
}

SparseMatrix<double>& DeformableObject::sparse_mass_matrix()
{
    return M_sparse;
}

SparseMatrix<double>& DeformableObject::sparse_inverse_mass_matrix()
{
    return M_inv_sparse;
}

MatrixXd& DeformableObject::mass_matrix()
{
    return M;
}

MatrixXd& DeformableObject::inverse_mass_matrix()
{
    return M_inv;
}


/////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////// PENALTY OBJECT /////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////

void PenaltyObject::wrapMesh(void)
{
    cMesh* mesh = this->getMesh(0);
    cTransform transform = mesh->getLocalTransform();

    int num_triangles = mesh->getNumTriangles();
    int num_vertices = mesh->getNumVertices();
    m_sensors.resize(num_triangles);
    vertexNeighbors.resize(num_vertices);
    vertexNormals.resize(num_vertices);

    for (int triangle_idx = 0; triangle_idx < num_triangles; triangle_idx++)
    {
         auto v1 = mesh->m_triangles->getVertexIndex0(triangle_idx);
         auto v2 = mesh->m_triangles->getVertexIndex1(triangle_idx);
         auto v3 = mesh->m_triangles->getVertexIndex2(triangle_idx);
         auto area = mesh->m_triangles->computeArea(triangle_idx);
         auto normal = mesh->m_triangles->computeNormal(triangle_idx, true);
         auto pos1 = mesh->m_vertices->getLocalPos(v1);
         auto pos2 = mesh->m_vertices->getLocalPos(v2);
         auto pos3 = mesh->m_vertices->getLocalPos(v3);
         auto localTransform = this->getLocalTransform();
         auto centroid = (pos1 + pos2 + pos3) / 3;
         int numpts = round(area * density);

        double r,g,b;
        generateRandomPastelColor(r,g,b);

        if (numpts == 1)
        {
            m_sensors[triangle_idx].emplace_back(new Sensor());

            // into the mesh
            cVector3d startPos = centroid - m_sensors[triangle_idx].back()->len*normal;
            cVector3d endPos = centroid ;//+ 0.75*m_sensors[triangle_idx].back()->len*normal;
            m_sensors[triangle_idx].back()->startPos = startPos.eigen();

            if (visualizeSensors)
            {
                m_rayViz.emplace_back(new cShapeLine(startPos, endPos));
                this->addChild(m_rayViz.back());

                m_rayViz.back()->m_colorPointA.set(r,g,b);
                m_rayViz.back()->m_colorPointB.set(r,g,b);
            }
        }

        if (numpts > 1)
        {
            auto pts_pos = distributePointsInTriangle(pos1, pos2, pos3, numpts);

            for (int r = 0; r < pts_pos.size(); r++)
            {
                m_sensors[triangle_idx].emplace_back(new Sensor());
                m_sensors[triangle_idx].back()->m_parent = this;

                // into the mesh
                cVector3d startPos = pts_pos[r] - 0.25*m_sensors[triangle_idx].back()->len*normal;
                cVector3d endPos = pts_pos[r] + 0.75*m_sensors[triangle_idx].back()->len*normal;

                m_sensors[triangle_idx].back()->startPos = startPos.eigen();
                m_sensors[triangle_idx].back()->normal = normal.eigen();

                if (visualizeSensors)
                {
                    m_rayViz.emplace_back(new cShapeLine(startPos, endPos));
                    this->m_parent->addChild(m_rayViz.back());
                    m_rayViz.back()->setLineWidth(2);
                    m_rayViz.back()->m_colorPointA.set(r,g,b);
                    m_rayViz.back()->m_colorPointB.set(r,g,b);
                }
            }
        }

    }

//    for (int i = 0 ; i < num_vertices; i++)
//    {
//        for (int triangle_idx = 0; triangle_idx < num_triangles; triangle_idx++) {
//            const auto &v1 = mesh->m_triangles->getVertexIndex0(triangle_idx);
//            const auto &v2 = mesh->m_triangles->getVertexIndex1(triangle_idx);
//            const auto &v3 = mesh->m_triangles->getVertexIndex2(triangle_idx);
//            if (i == v1 || i == v2 || i == v3)
//            {
//                if (i != v1)
//                    vertexNeighbors[i].insert(v1);
//                if (i != v2)
//                    vertexNeighbors[i].insert(v2);
//                if (i != v3)
//                    vertexNeighbors[i].insert(v3);
//            }
//        }
//    }
//
//    for (int i = 0 ; i < num_vertices; i++)
//    {
//        Vector3d normalAvg(0,0,0);
//        int div = 0;
//        for (auto triangle_idx : vertexNeighbors[i])
//        {
//            auto v1 = mesh->m_triangles->getVertexIndex0(triangle_idx);
//            auto v2 = mesh->m_triangles->getVertexIndex1(triangle_idx);
//            auto v3 = mesh->m_triangles->getVertexIndex2(triangle_idx);
//            auto normal = mesh->m_triangles->computeNormal(triangle_idx, true);
//            normalAvg += normal.eigen();
//            div++;
//        }
//        normalAvg /= div;
//    }

}

void PenaltyObject::updateSensors()
{
    Affine3d transform;
    transform.translation() = this->x;
    transform.linear() = this->q.toRotationMatrix();
    int i = 0;

    for (auto triangle_sensors : m_sensors)
    {
        for (auto sensor : triangle_sensors)
        {
            sensor->globalStartPos = transform*sensor->startPos;
            sensor->globalNormal = transform*sensor->normal;
            m_rayViz[i]->setLocalPos(sensor->globalStartPos);
            i++;
        }
    }
}