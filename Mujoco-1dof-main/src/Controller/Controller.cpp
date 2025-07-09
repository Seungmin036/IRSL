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
    rotor_inertia_matrix_.resize(nq,nq); rotor_inertia_matrix_.setZero(); rotor_inertia_matrix_.diagonal() << 0.4, 0.4, 0.4, 0.4, 0.2, 0.2, 0.2   ;

    
    // LuGre parameter initialization
    z     = Eigen::VectorXd::Zero(nq);  
    z_prev     = Eigen::VectorXd::Zero(nq);  
    sig_0 = Eigen::VectorXd::Constant(nq, 275.0);
    sig_1 = Eigen::VectorXd::Constant(nq, 45.2);
    sig_2 = Eigen::VectorXd::Constant(nq, 1.819);
    Fc    = Eigen::VectorXd::Constant(nq, 6.975);
    Fs    = Eigen::VectorXd::Constant(nq, 8.875);
    vs    = Eigen::VectorXd::Constant(nq, 0.06109);
    

    theta_n_prev = Eigen::VectorXd::Zero(nq);
    dtheta_n_prev = Eigen::VectorXd::Zero(nq);
    sigma_prev = Eigen::VectorXd::Zero(nq);
    tau_f_prev = Eigen::VectorXd::Zero(nq);
    tau_j_prev = Eigen::VectorXd::Zero(nq);
    tau_c_prev = Eigen::VectorXd::Zero(nq);
    
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
        // u = u - lugre_friction(robot_state);
        break;
    case CONTROLLER_SELECTOR::JOINT_PD_FRIC:
        u=this->PD_controller_AS_GC(robot_state);
        // u = u - lugre_friction(robot_state);
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
    //Kp.diagonal() << 1000, 1000, 1000, 1000, 500, 500, 500;
    //Kd.diagonal() << 250, 250, 250, 250, 150, 150, 150; 
    Kp.diagonal() << 200, 200, 200, 200, 100, 100, 100;
    Kd.diagonal() << 50, 50, 50, 50, 30, 30, 30; 

    Eigen::VectorXd control_motor_torque;
    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque.setZero();

    gravity_compensation_torque = robot.GetGravity(robot_state.q);

    theta_des = robot_state.q_d+joint_stiffness_matrix_.inverse()*gravity_compensation_torque;
    control_motor_torque = -Kp*(robot_state.theta-theta_des)-Kd*(robot_state.dtheta - robot_state.dq_d)+gravity_compensation_torque;

    return control_motor_torque;
}

Eigen::VectorXd Controller::lugre_friction(const RobotState & robot_state)
{
    Eigen::VectorXd friction(nq);
 

    for (int i = 0; i < nq; ++i) {
        double g = Fc(i) + (Fs(i) - Fc(i)) * std::exp(-std::pow(robot_state.dtheta(i) / vs(i), 2));
        // double dz_i = robot_state.dtheta(i) - (sig_0(i) * std::abs(robot_state.dtheta(i)) / g) * z_prev(i);
        // z(i) = z_prev(i) + dz_i * step_time_;
        double dz_i = robot_state.dtheta(i) - (sig_0(i) * std::abs(robot_state.dtheta(i)) / g) * z(i);
        z(i) = z(i) + dz_i * step_time_;


        // calculate friction
        friction(i) = sig_0(i) * z(i) + sig_1(i) * dz_i + sig_2(i) * robot_state.dtheta(i);
        
    }
    // std::cout << "tau_f = " << friction.transpose() << std::endl;
    return friction;
}

Eigen::VectorXd Controller::friction_observer_PD(const RobotState & robot_state)  
{
    // Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( tau_c_prev - tau_j_prev + sigma_prev - tau_f_prev);
    // Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( tau_c_prev - tau_j_prev);

    Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( tau_c_prev - robot_state.tau_J);
    // std::cout << "ddtheta_n" << ddtheta_n.transpose() << std::endl;

    Eigen::VectorXd dtheta_n = dtheta_n_prev + step_time_ * ddtheta_n;
    std::cout << "dtheta_n" << dtheta_n.transpose() << std::endl;
    Eigen::VectorXd theta_n = theta_n_prev + step_time_ * dtheta_n_prev;    
    // std::cout << "[DEBUG] tau_c_prev = " << tau_c_prev.transpose() << std::endl;
    // std::cout << "[DEBUG] tau_J      = " << robot_state.tau_J.transpose() << std::endl;
    // std::cout << "[DEBUG] ddtheta_n  = " << ddtheta_n.transpose() << std::endl;

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
    // update
    theta_n_prev = theta_n;
    dtheta_n_prev = dtheta_n;
    sigma_prev = sigma;
    tau_f_prev = tau_f;
    tau_j_prev = robot_state.tau_J;
    // std::cout << "tau_f_hat = " << tau_f.transpose() << std::endl;

    return tau_f;
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
    Eigen::VectorXd tau_f_hat = friction_observer_PD(robot_state);

    //update tau_c
    tau_c_prev = tau_c;

    // Eigen::VectorXd control_motor_torque = tau_c-tau_f_hat;
    Eigen::VectorXd control_motor_torque = tau_c;

    return control_motor_torque;
}