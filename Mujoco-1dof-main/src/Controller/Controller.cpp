#include <Controller/Controller.hpp>

Controller::Controller(const pinocchio::Model & pinocchio_model) : robot(pinocchio_model)
{
    pinocchio_model_ = pinocchio_model;
    nq = pinocchio_model_.nq;
    nv = pinocchio_model_.nv;
}   

void Controller::InitController(const RobotState & robot_state_init)
{
    joint_stiffness_matrix_.resize(nq,nq); joint_stiffness_matrix_.setZero();
    rotor_inertia_matrix_.resize(nq,nq); rotor_inertia_matrix_.setZero();

    //for initial model parameter

    // LuGre parameter initialization
    z     = Eigen::VectorXd::Zero(nq);
    sig_0 = Eigen::VectorXd::Constant(nq, 2750.0);
    sig_1 = Eigen::VectorXd::Constant(nq, 45.2);
    sig_2 = Eigen::VectorXd::Constant(nq, 1.819);
    Fc    = Eigen::VectorXd::Constant(nq, 6.975);
    Fs    = Eigen::VectorXd::Constant(nq, 8.875);
    vs    = Eigen::VectorXd::Constant(nq, 0.06109);

}

Eigen::VectorXd Controller::GetControlInput(const RobotState & robot_state, const int& selector)
{   
    Eigen::VectorXd u;
    switch (selector)
    {
    case CONTROLLER_SELECTOR::JOINT_PD:
        u=this->PD_controller(robot_state);
        break;
    case CONTROLLER_SELECTOR::JOINT_PD_FRIC:
        // u=this->PD_controller_AS_GC(robot_state);
        break;

    case CONTROLLER_SELECTOR::TASK_PD:
        // u=this->PD_controller_AS_GC_task_space(robot_state);
        break;

    case CONTROLLER_SELECTOR::TASK_PD_FRIC:
        // u=this->PD_controller_AS_GC_task_space(robot_state);
        break;
    
    default:
        break;
    }

    return u;
}

Eigen::VectorXd Controller::PD_controller(const RobotState & robot_state)
{
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity();
    Kd.setIdentity();

    // Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kp = Kp*50;
    // Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; //PD controller derivative gain
    Kd = Kd*3;

    Eigen::VectorXd control_motor_torque;

    if(robot_state.is_rigid)
    {   
        Eigen::VectorXd gravity_compensation_torque(nv);
        gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
        control_motor_torque = -Kp*(robot_state.q-robot_state.q_d)-Kd*(robot_state.dq - robot_state.dq_d)+gravity_compensation_torque;
    }
    else
    {
        Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
        gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
        // std::cout<<"gravity_compensation_torque : \n"<<gravity_compensation_torque.transpose()<<std::endl;

        theta_des = robot_state.q_d+joint_stiffness_matrix_.inverse()*gravity_compensation_torque;
        control_motor_torque = -Kp*(robot_state.theta-theta_des)-Kd*(robot_state.dtheta - robot_state.dq_d)+gravity_compensation_torque;
        control_motor_torque -=lugre_fricion(const Robotstate & robot_state)
    }
    
    return control_motor_torque;
}

Eigen::VectorXd Controller::lugre_friction(const RobotState & robot_state)
{
    Eigen::VectorXd friction(nq);  

    for (int i = 0; i < nq; ++i) {
        double v = robot_state.dtheta(i);
        double g = Fc(i) + (Fs(i) - Fc(i)) * std::exp(-std::pow(v / vs(i), 2));
        double dz_i = v - (sig_0(i) * std::abs(v) / g) * z(i);
        z(i) += dz_i * step_time_;

        // calculate friction
        friction(i) = sig_0(i) * z(i) + sig_1(i) * dz_i + sig_2(i) * v;
    }

    return friction;
}