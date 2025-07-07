#include <Controller/Controller.hpp>

Controller::Controller(const pinocchio::Model & pinocchio_model) : robot(pinocchio_model)
{
    pinocchio_model_ = pinocchio_model;
    nq = pinocchio_model_.nq;
    nv = pinocchio_model_.nv;
}   

void Controller::InitController(const RobotState & robot_state_init)
{
    //for initial model parameter
    joint_stiffness_matrix_.resize(nq,nq); joint_stiffness_matrix_.setZero(); joint_stiffness_matrix_.diagonal() << 8000, 8000, 8000, 8000, 7000, 7000, 7000;
    rotor_inertia_matrix_.resize(nq,nq); rotor_inertia_matrix_.setZero(); rotor_inertia_matrix_.diagonal() << 150,150,150,150,120,120,120   ;

    
    // LuGre parameter initialization
    z     = Eigen::VectorXd::Zero(nq);
    sig_0 = Eigen::VectorXd::Constant(nq, 2750.0);
    sig_1 = Eigen::VectorXd::Constant(nq, 45.2);
    sig_2 = Eigen::VectorXd::Constant(nq, 1.819);
    Fc    = Eigen::VectorXd::Constant(nq, 6.975);
    Fs    = Eigen::VectorXd::Constant(nq, 8.875);
    vs    = Eigen::VectorXd::Constant(nq, 0.06109);

    // observer gain
    tau_f_hat_ = Eigen::VectorXd::Zero(nq);
    e_nr_ = Eigen::VectorXd::Zero(nq);
    e_dot_nr_ = Eigen::VectorXd::Zero(nq);
    theta_nom_ = robot_state_init.theta;
    theta_dot_nom_ = robot_state_init.dtheta;

    Gamma_ = Eigen::MatrixXd::Identity(nq, nq) * 10000;
    Gamma_p_ = Eigen::MatrixXd::Identity(nq, nq) * 100000*step_time_;
    K_lpf_= Eigen::VectorXd::Zero(nq);
    K_lpf_ << 150, 150, 150, 150, 120, 120, 120;
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
        u=this->PD_controller_AS_GC(robot_state);
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

Eigen::VectorXd Controller::quasi_static_estimate_q(const RobotState & robot_state)
{
    Eigen::VectorXd gravity = robot.GetGravity(robot_state.q_d);
    Eigen::VectorXd q_est = robot_state.theta - joint_stiffness_matrix_.inverse() * gravity;

    return q_est;
}

Eigen::VectorXd Controller::PD_controller(const RobotState & robot_state)
{
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity();  Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kd.setIdentity();  Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; 

    Eigen::VectorXd control_motor_torque;
    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque = robot.GetGravity(robot_state.q_d);

    theta_des = robot_state.q_d+joint_stiffness_matrix_.inverse()*gravity_compensation_torque;
    control_motor_torque = -Kp*(robot_state.theta-theta_des)-Kd*(robot_state.dtheta - robot_state.dq_d)+gravity_compensation_torque;
    
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

Eigen::VectorXd Controller::friction_observer_PD(const RobotState & robot_state, Eigen::VectorXd control_motor_torque)
{
    Eigen::VectorXd theta_ddot_nom = rotor_inertia_matrix_.inverse() * (control_motor_torque - robot_state.tau_J);

    theta_dot_nom_ += step_time_ * theta_ddot_nom;
    theta_nom_     += step_time_ * theta_dot_nom_;

    e_nr_ = theta_nom_ - robot_state.theta;
    e_dot_nr_ = theta_dot_nom_ - robot_state.dtheta;

    Eigen::VectorXd sigma_hat = -rotor_inertia_matrix_ * Gamma_ * (e_dot_nr_ + Gamma_p_ * e_nr_);
    for (int i = 0; i < nq; ++i)
    {
        double alpha = std::exp(-K_lpf_(i) * step_time_);
        tau_f_hat_(i) = alpha * tau_f_hat_(i) + (1.0 - alpha) * sigma_hat(i);
    }

    return tau_f_hat_;
}


Eigen::VectorXd Controller::PD_controller_AS_GC(const RobotState & robot_state)
{
    // PD controller with active stiffness + gravity compensation + friction observer
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity();
    Kd.setIdentity();

    Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; 
    

    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
    theta_des = robot_state.q_d + joint_stiffness_matrix_.inverse() * gravity_compensation_torque;

    Eigen::VectorXd tau_c = -Kp * (robot_state.theta - theta_des)- Kd * (robot_state.dtheta - robot_state.dq_d)+ gravity_compensation_torque;

    Eigen::VectorXd tau_f_hat = friction_observer_PD(robot_state, tau_c);
    Eigen::VectorXd control_motor_torque = tau_c - tau_f_hat;
    control_motor_torque += lugre_friction(robot_state);

    return control_motor_torque;
}