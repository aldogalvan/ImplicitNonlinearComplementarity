
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP

#include <Eigen/Dense>

using namespace Eigen;

inline double computeTetrahedronVolume(const Vector3d& vertex1, const Vector3d& vertex2, const Vector3d& vertex3, const Vector3d& vertex4)
{
    // Calculate vectors representing three edges of the tetrahedron
    Vector3d edge1 = vertex2 - vertex1;
    Vector3d edge2 = vertex3 - vertex1;
    Vector3d edge3 = vertex4 - vertex1;

    // Calculate the volume using the scalar triple product
    double volume = edge1.dot(edge2.cross(edge3)) / 6.0;

    return std::abs(volume);
}

inline Eigen::Quaterniond generateQuaternion(double rotation_x_deg, double rotation_y_deg, double rotation_z_deg)
{
    double rotation_x_rad = rotation_x_deg * M_PI / 180.0; // Convert degrees to radians
    double rotation_y_rad = rotation_y_deg * M_PI / 180.0;
    double rotation_z_rad = rotation_z_deg * M_PI / 180.0;

    Eigen::AngleAxisd rot_x(rotation_x_rad, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd rot_y(rotation_y_rad, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd rot_z(rotation_z_rad, Eigen::Vector3d::UnitZ());

    Eigen::Quaterniond quaternion = rot_z * rot_y * rot_x;

    return quaternion;
}

inline Eigen::Quaterniond multiplyQuaternionScalar( double scalar,const Eigen::Quaterniond& q)
{
    Eigen::Quaterniond result;
    result.coeffs() = q.coeffs() * scalar;
    return result;
}


inline Eigen::Quaterniond sumQuaternions(const Eigen::Quaterniond& q1, const Eigen::Quaterniond& q2)
{
    Eigen::Quaterniond sum;
    sum.coeffs() = q1.coeffs() + q2.coeffs();
    return sum;
}


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
