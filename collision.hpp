
#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP

#include <Eigen/Dense>
#include <vector>
#include "aabb.hpp"

using namespace Eigen;
using namespace std;

struct ColInfo {
    // time of collision
    double t;
    // collision normal
    Vector3d normal;
    // collision depth
    double depth;
    // collision body indices
    int bodyIdxA, bodyIdxB;
};

bool findCollisions(MatrixXd object1_start, MatrixXd object1_end, MatrixXi object1_tris,
                    MatrixXd object2_start, MatrixXd object2_end, MatrixXi object2_tris,
                    vector<ColInfo*>& collisions);

#endif //IMPLICITNONLINEARCOMPLEMENTARITY_COLLISION_HPP
