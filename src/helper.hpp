
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP

#include <Eigen/Dense>

using namespace Eigen;

inline Vector3d checkNearZero(const Eigen::Vector3d& vector) {
    double threshold = 1e-6;

    if (vector.norm() < threshold) {
        return Eigen::Vector3d::Zero();
    } else {
        return vector;
    }
}

inline void cutoff_negative_values(Eigen::VectorXd& v) {
    for (int i = 0; i < v.size(); i++) {
        if (v(i) < 0) {
            v(i) = 0;
        }
    }
}

inline Eigen::MatrixXd vec2skew(const Eigen::Vector3d& vec)
{
    Eigen::MatrixXd skew(3, 3);
    skew <<  0, -vec.z(), vec.y(),
            vec.z(), 0, -vec.x(),
            -vec.y(), vec.x(), 0;
    return skew;
}

inline void applyQuaternionRotation(Eigen::MatrixXd& vertices, const Eigen::Quaterniond& rotation)
{
    // Get the number of vertices and the dimensionality of each vertex
    const int numVertices = vertices.rows();
    const int vertexDim = vertices.cols();

    // Apply the rotation to each vertex
    for (int i = 0; i < numVertices; ++i)
    {
        // Extract the vertex coordinates
        Eigen::Vector3d vertex = vertices.row(i).head<3>();

        // Apply the rotation to the vertex coordinates
        vertex = rotation * vertex;

        // Update the vertex in the vertices matrix
        vertices.row(i).head<3>() = vertex;
    }
}

static Matrix3d anglesToRotationMatrix(Vector3d angles)
{

    double angleX = angles(0); double angleY = angles(1); double angleZ = angles(2);
    Matrix3d rotationMatrix;

    rotationMatrix = AngleAxisd(angleX, Vector3d::UnitX())
                     * AngleAxisd(angleY, Vector3d::UnitY())
                     * AngleAxisd(angleZ, Vector3d::UnitZ());

    return rotationMatrix;
}

static Vector3d rotationMatrixToAngles(const Matrix3d& rotationMatrix)
{
    Vector3d euler = rotationMatrix.eulerAngles(0, 1, 2); // Angles in XYZ order

    return Vector3d(euler[0],euler[1],euler[2]);
}

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP
