#pragma once

#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/fwd.hpp>
#include "pinocchio/algorithm/joint-configuration.hpp"

#include "sensor_msgs/JointState.h"
#include "geometry_msgs/PoseArray.h"
#include "std_msgs/String.h"

#include "ros/ros.h"
#include "vector"
#include "mjxmacro.h"
#include "uitools.h"
#include <Eigen/Dense>
#include <data/RobotState.hpp>
#include <mujoco/mujoco.h>

#include "data/ObjectState.hpp"
#include <Controller/Controller.hpp>


class ControllerNode
{
private:
    ros::NodeHandle n;

    const int nq;
    const int nv;
    const int nu;

    ros::Subscriber object_name_msg_sub;

    //For publish current joint state
    const std::string frame_id="world";
    std::vector<std::string> joint_names;
    ros::Publisher joint_state_msg_pub;
    ros::Publisher object_state_msg_pub;
    std::vector<std::string> object_names;

    //controller
    Controller controller;

    //Member varialbe
    RobotState robot_state;
    RobotState robot_state_init;
    ObjectState object_state;
    const bool is_rigid;


public:
    ControllerNode(ros::NodeHandle &nh, const pinocchio::Model &pinnochio_model);
    ~ControllerNode();

    void InitControllerNode(const mjModel* m, mjData* d);
    void UpdateRobotState(const mjModel* m, mjData* d, RobotState & robot_state);
    void Control_Loop(const mjModel* m, mjData* d);
};