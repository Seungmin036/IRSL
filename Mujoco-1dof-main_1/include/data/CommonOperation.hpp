#pragma once

#include <Eigen/Core>

inline Eigen::Matrix3d skew(Eigen::Vector3d v)
{
    Eigen::Matrix3d out;

    out << 0,-v(2),v(1),
            v(2),0, -v(0),
            -v(1),v(0),0;

    return out;
}

inline Eigen::Matrix3d RPY2RotMatrix(const double &roll, const double &pitch, const double &yaw)
{
    Eigen::AngleAxisd rollAngle(roll, Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd pitchAngle(pitch, Eigen::Vector3d::UnitY());
    Eigen::AngleAxisd yawAngle(yaw, Eigen::Vector3d::UnitZ());

    Eigen::Quaternion<double> q =  yawAngle *pitchAngle * rollAngle;

    return q.matrix();
}