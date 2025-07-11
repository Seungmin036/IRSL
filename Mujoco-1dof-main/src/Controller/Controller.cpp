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
    rotor_inertia_matrix_.resize(nq,nq); rotor_inertia_matrix_.setZero(); rotor_inertia_matrix_.diagonal() << 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2 ;

    
    // LuGre parameter initialization
    z     = Eigen::VectorXd::Zero(nq);  
    sig_0 = Eigen::VectorXd::Constant(nq, 2750.0);
    sig_1 = Eigen::VectorXd::Constant(nq, 45.2);
    sig_2 = Eigen::VectorXd::Constant(nq, 1.819);
    Fc    = Eigen::VectorXd::Constant(nq, 6.975);
    Fs    = Eigen::VectorXd::Constant(nq, 8.875);
    vs    = Eigen::VectorXd::Constant(nq, 0.06109);
    

    theta_n_prev = robot_state_init.theta;
    dtheta_n_prev = robot_state_init.dtheta;
    sigma_prev = Eigen::VectorXd::Zero(nq);
    tau_f_prev = Eigen::VectorXd::Zero(nq);
    tau_j_prev = Eigen::VectorXd::Zero(nq);
    u_prev = Eigen::VectorXd::Zero(nq);
    
    // observer gain
    Gamma_ = Eigen::MatrixXd::Identity(nq, nq) * 50000;
    Gamma_p_ = 10*step_time_*Gamma_;
    K_lpf_= Eigen::VectorXd::Zero(nq); K_lpf_ << 150, 150, 150, 150, 120, 120, 120;
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
        // u=this->PD_controller_AS_GC(robot_state)+lugre_friction(robot_state) - friction_observer_PD(robot_state, u_prev);
        u=this->PD_controller_AS_GC(robot_state) - lugre_friction(robot_state);
        u_prev = u;
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
    Kp.setIdentity();  Kd.setIdentity(); 
    Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; 

    Eigen::VectorXd control_motor_torque;
    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);

    gravity_compensation_torque = robot.GetGravity(robot_state.q);

    theta_des = robot_state.q_d+joint_stiffness_matrix_.inverse()*gravity_compensation_torque;
    control_motor_torque = -Kp*(robot_state.theta-theta_des)-Kd*(robot_state.dtheta - robot_state.dq_d)+gravity_compensation_torque;

    return control_motor_torque;
}

Eigen::VectorXd Controller::lugre_friction(const RobotState & robot_state)
{
    Eigen::VectorXd friction(nq);
    for (int i = 0; i < nq; ++i) {
        double dtheta_i = robot_state.dtheta(i);
        // 작은 속도 threshold 적용 
        if (std::abs(dtheta_i) < 1e-5)
            dtheta_i = 0.0;

        double g = Fc(i) + (Fs(i) - Fc(i)) * std::exp(-std::pow(dtheta_i / vs(i), 2));
        double dz_i = dtheta_i - sig_0(i) * (std::abs(dtheta_i) / g) * z(i);
        z(i) += dz_i * step_time_;
        friction(i) = sig_0(i) * z(i) + sig_1(i) * dz_i + sig_2(i) * dtheta_i;
    }
    return friction;
}

Eigen::VectorXd Controller::friction_observer_PD(const RobotState & robot_state, VectorXd Control_input)  
{
    Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( Control_input - tau_j_prev + sigma_prev - tau_f_prev );   
    Eigen::VectorXd dtheta_n = dtheta_n_prev + step_time_ * ddtheta_n;
    Eigen::VectorXd theta_n = theta_n_prev + step_time_ * dtheta_n; 
    Eigen::VectorXd de_nr = dtheta_n - robot_state.dtheta;
    Eigen::VectorXd e_nr = theta_n - robot_state.theta;

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nq, nq);
    Eigen::VectorXd dde_nr = -(I + Gamma_ * step_time_).inverse() * Gamma_ * ((I + Gamma_p_ * step_time_) * de_nr + Gamma_p_ * e_nr);
    Eigen::VectorXd sigma = -rotor_inertia_matrix_ * Gamma_ * (step_time_ * dde_nr + (I + Gamma_p_ * step_time_) * de_nr + Gamma_p_ * e_nr);
    Eigen::VectorXd tau_f(nq);
    for (int i = 0; i < nq; ++i) {
       double alpha_i = std::exp(-K_lpf_(i) * step_time_);
       tau_f(i) = -alpha_i * tau_f_prev(i) + (1.0 - alpha_i) * sigma(i);
    }

    theta_n_prev = theta_n;
    dtheta_n_prev = dtheta_n;
    sigma_prev = sigma;
    tau_f_prev = tau_f;
    tau_j_prev = robot_state.tau_J;
    return tau_f;
}

Eigen::VectorXd Controller::PD_controller_AS_GC(const RobotState & robot_state)
{
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity();
    Kd.setIdentity();

    Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; 
    

    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
    theta_des = robot_state.q_d + joint_stiffness_matrix_.inverse() * gravity_compensation_torque;

    Eigen::VectorXd tau_c = -Kp * (robot_state.theta - theta_des)- Kd * (robot_state.dtheta - robot_state.dq_d)+ gravity_compensation_torque;
    Eigen::VectorXd control_motor_torque = tau_c;

    return control_motor_torque;
}