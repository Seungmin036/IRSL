#pragma once

#include <Eigen/Core>

struct ObjectState
{
    int num_object = 0;
    
    std::vector<Eigen::Vector3d> positions;
    std::vector<Eigen::Matrix3d> orientations;
    std::vector<Eigen::Quaterniond> quat;

};