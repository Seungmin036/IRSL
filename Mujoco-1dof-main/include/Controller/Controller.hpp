#pragma once

#pragma once
#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/fwd.hpp>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <chrono>
#include <thread>

#include "ros/ros.h"

#include "Dynamics/Dynamics.hpp"
#include "data/RobotState.hpp"
#include <data/CommonOperation.hpp>

using namespace std;
using namespace Eigen;

class Controller
{
    
    private:
    int nq;
    int nv;

    double step_time_=0.001; // 1000Hz defined in franka_panda__flexible.xml file.

    // for LuGre 
    Eigen::VectorXd z, z_prev;
    Eigen::VectorXd sig_0, sig_1, sig_2;
    Eigen::VectorXd Fc, Fs, vs;
    

    //for flexible joint robot parameter
    Eigen::MatrixXd joint_stiffness_matrix_;
    Eigen::MatrixXd rotor_inertia_matrix_;

    //for observer
    Eigen::VectorXd theta_n;
    Eigen::VectorXd dtheta_n;
    Eigen::VectorXd sigma_prev;
    Eigen::VectorXd tau_f_prev;
    Eigen::VectorXd u_prev;
    Eigen::VectorXd tau_j_prev;

    Eigen::MatrixXd Gamma_, Gamma_p_;     
    Eigen::VectorXd K_lpf_;
    Eigen::MatrixXd L, Lp;  

    //for RNEA 
    Dynamics::articulated_system robot;
    pinocchio::Model pinocchio_model_;

    ros::NodeHandle nh;
    
    //*****************************Controller**********************************//

    Eigen::VectorXd quasi_static_estimate_q(const RobotState & robot_state);

    Eigen::VectorXd lugre_friction(const RobotState & robot_state);

    Eigen::VectorXd friction_observer_PD(const RobotState & robot_state, VectorXd& Control_input);
    
    Eigen::VectorXd friction_observer_L1_PD(const RobotState & robot_state, VectorXd& Control_input);

    Eigen::VectorXd PD_controller(const RobotState & robot_state);

    Eigen::VectorXd PD_controller_AS_GC(const RobotState & robot_state);

    Eigen::VectorXd PD_controller_AS_GC_task_space(const RobotState & robot_state);

    // Eigen::MatrixXd Null_space_projection(const RobotState & robot_state);

    //*************************************************************************//

    public:
    int loop_index_=0;

    Controller(const pinocchio::Model & pinocchio_model);

    void InitController(const RobotState & robot_state_init);
    Eigen::VectorXd GetControlInput(const RobotState & robot_state, const int& selectors);

};