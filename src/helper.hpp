
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_HELPER_HPP

#include <Eigen/Dense>

using namespace Eigen;

//! A bunch of useful functions

inline Quaterniond angularVelocityToQuaternion(const Vector3d& angular_velocity)
{
    Quaterniond angular_velocity_quaternion;
    angular_velocity_quaternion.w() = 0.0;
    angular_velocity_quaternion.vec() = angular_velocity;

    return angular_velocity_quaternion;
}

inline bool is_ill_conditioned(const MatrixXd& matrix, double threshold) {
    JacobiSVD<MatrixXd> svd(matrix);
    double condNumber = svd.singularValues()(0) / svd.singularValues()(svd.singularValues().size() - 1);
    return condNumber > threshold;
}


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

inline Quaterniond generateQuaternion(double rotation_x_deg, double rotation_y_deg, double rotation_z_deg)
{
    double rotation_x_rad = rotation_x_deg * M_PI / 180.0; // Convert degrees to radians
    double rotation_y_rad = rotation_y_deg * M_PI / 180.0;
    double rotation_z_rad = rotation_z_deg * M_PI / 180.0;

    AngleAxisd rot_x(rotation_x_rad, Vector3d::UnitX());
    AngleAxisd rot_y(rotation_y_rad, Vector3d::UnitY());
    AngleAxisd rot_z(rotation_z_rad, Vector3d::UnitZ());

    Quaterniond quaternion = rot_z * rot_y * rot_x;

    return quaternion;
}

inline Quaterniond multiplyQuaternionScalar( double scalar,const Quaterniond& q)
{
    Quaterniond result;
    result.coeffs() = q.coeffs() * scalar;
    return result;
}


inline Quaterniond sumQuaternions(const Quaterniond& q1, const Quaterniond& q2)
{
    Quaterniond sum;
    sum.coeffs() = q1.coeffs() + q2.coeffs();
    return sum;
}


inline Vector3d checkNearZero(const Vector3d& vector) {
    double threshold = 1e-6;

    if (vector.norm() < threshold) {
        return Vector3d::Zero();
    } else {
        return vector;
    }
}

inline void cutoff_negative_values(VectorXd& v) {
    for (int i = 0; i < v.size(); i++) {
        if (v(i) < 0) {
            v(i) = 0;
        }
    }
}

inline MatrixXd vec2skew(const Vector3d& vec)
{
    MatrixXd skew(3, 3);
    skew <<  0, -vec.z(), vec.y(),
            vec.z(), 0, -vec.x(),
            -vec.y(), vec.x(), 0;
    return skew;
}

inline void applyQuaternionRotation(MatrixXd& vertices, const Quaterniond& rotation)
{
    // Get the number of vertices and the dimensionality of each vertex
    const int numVertices = vertices.rows();
    const int vertexDim = vertices.cols();

    // Apply the rotation to each vertex
    for (int i = 0; i < numVertices; ++i)
    {
        // Extract the vertex coordinates
        Vector3d vertex = vertices.row(i).head<3>();

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
