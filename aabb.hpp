

#ifndef IMPLICITNONLINEARCOMPLEMENTARITY_AABB_HPP
#define IMPLICITNONLINEARCOMPLEMENTARITY_AABB_HPP

#include <Eigen/Core>
#include <memory>
#include <vector>

struct AABBNode;
struct BBox;

std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> startPos,
                                    Eigen::Ref<const Eigen::MatrixXd> endPos,
                                    Eigen::Ref<const Eigen::MatrixXi> F);

std::shared_ptr<AABBNode> buildAABB(Eigen::Ref<const Eigen::MatrixXd> pos,
                                    Eigen::Ref<const Eigen::MatrixXi> F );

struct Collision {
    int collidingTriangle1; // Triangle from the "left" AABBNode
    int collidingTriangle2; // Triangle from the "right" AABBNode
};

// Intersect AABB "left" vs "right"
void intersect(std::shared_ptr<AABBNode> left,
               std::shared_ptr<AABBNode> right,
               std::vector<Collision>& collisions);


#endif //IMPLICITNONLINEARCOMPLEMENTARITY_AABB_HPP
