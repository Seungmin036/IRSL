#pragma once

#include <Eigen/Core>

struct RobotState
{
  Eigen::VectorXd q; // link position
  Eigen::VectorXd dq; // qdot

  Eigen::VectorXd q_d; // desired link positions
  Eigen::VectorXd dq_d; // desired link velocity

  Eigen::VectorXd theta; // motor position
  Eigen::VectorXd dtheta;

  Eigen::VectorXd tau_J; // JTS measurement

  std::vector<Eigen::Affine3d> ee_htm; // End-effector position
  std::vector<Eigen::Affine3d> ee_htm_d;

  bool is_rigid;

  void print()
  { 
    std::cout<<"----------------------------------------"<<std::endl;
    std::cout<<"q : " << q.transpose()<<std::endl;
    std::cout<<"dq : " << dq.transpose()<<std::endl;
    std::cout<<"theta : " << theta.transpose()<<std::endl;
    std::cout<<"dtheta : " << dtheta.transpose()<<std::endl;
    std::cout<<"q_d : " << q_d.transpose()<<std::endl;
    std::cout<<"dq_d : " << dq_d.transpose()<<std::endl;

  }
};

enum CONTROLLER_SELECTOR
{
  JOINT_PD,
  JOINT_PD_FRIC,
  TASK_PD,
  TASK_PD_FRIC,
};