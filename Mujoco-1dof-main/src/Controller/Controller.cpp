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
    

    theta_n = robot_state_init.theta;
    dtheta_n = robot_state_init.dtheta;
    sigma_prev = Eigen::VectorXd::Zero(nq);
    tau_f_prev = Eigen::VectorXd::Zero(nq);
    tau_j_prev = Eigen::VectorXd::Zero(nq);
    u_prev = Eigen::VectorXd::Zero(nq);
    
    // observer gain
    Gamma_ = Eigen::MatrixXd::Identity(nq, nq) * 1e+4;
    Gamma_p_ = 100*Eigen::MatrixXd::Identity(7,7);
    K_lpf_= Eigen::VectorXd::Zero(nq); K_lpf_ << 150, 150, 150, 150, 120, 120, 120;
    L = 80*Eigen::MatrixXd::Identity(nq, nq);
    Lp = 40*Eigen::MatrixXd::Identity(7,7);
}

Eigen::VectorXd Controller::GetControlInput(const RobotState & robot_state, const int& selector)
{   
    Eigen::VectorXd u;
    switch (selector)
    {
    case CONTROLLER_SELECTOR::JOINT_PD:
    {   Eigen::VectorXd lugre = lugre_friction(robot_state);
        Eigen::VectorXd tau_f_obs = friction_observer_PD(robot_state, u_prev);
        u_prev = PD_controller(robot_state);
        u= u_prev - lugre - tau_f_obs;
        std::ofstream log("friction_compare_log.csv", std::ios::app);

        // 루그레 마찰력
        for (int i = 0; i < nq; ++i) log << lugre(i) << ",";

        // 옵서버 마찰 추정값
        for (int i = 0; i < nq; ++i) {
            log << tau_f_obs(i);
            if (i != nq - 1) log << ",";
            else log << "\n";
        }

        log.close();
        break;
    }
    case CONTROLLER_SELECTOR::JOINT_PD_FRIC:
    {   Eigen::VectorXd lugre = lugre_friction(robot_state);
        Eigen::VectorXd tau_f_obs = friction_observer_L1_PD(robot_state, u_prev);
        u_prev = PD_controller_AS_GC(robot_state);
        u= u_prev - lugre - tau_f_obs;
        std::ofstream log("friction_compare_log.csv", std::ios::app);

        // 루그레 마찰력
        for (int i = 0; i < nq; ++i) log << lugre(i) << ",";

        // 옵서버 마찰 추정값
        for (int i = 0; i < nq; ++i) {
            log << tau_f_obs(i);
            if (i != nq - 1) log << ",";
            else log << "\n";
        }

        log.close();
        break;
    }

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
    Kp.diagonal() << 50,50,50,50,50,50,50;
    Kd.diagonal() << 3,3,3,3,3,3,3;

    Eigen::VectorXd control_motor_torque;
    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);

    gravity_compensation_torque = robot.GetGravity(robot_state.q);

    theta_des = robot_state.q_d+joint_stiffness_matrix_.inverse()*gravity_compensation_torque;
    control_motor_torque = -Kp*(theta_n-theta_des)-Kd*(dtheta_n - robot_state.dq_d)+gravity_compensation_torque;

    return control_motor_torque;
}

Eigen::VectorXd Controller::lugre_friction(const RobotState & robot_state)
{
    Eigen::VectorXd friction(nq);
    Eigen::VectorXd dz(nq);

    for (int i = 0; i < nq; ++i)
    {
        double dtheta = robot_state.dtheta(i);
        double g = Fc(i) + (Fs(i) - Fc(i)) * std::exp(-std::pow(dtheta / vs(i), 2));

        // 안정화: g의 최소값 제한
        g = std::max(g, 1e-3);

        if (std::abs(dtheta) > 1e-6)
            dz(i) = dtheta - sig_0(i)*(std::abs(dtheta) / g) * z(i);
        else
            dz(i) = 0.0;

        z(i) += dz(i) * step_time_;

        friction(i) = sig_0(i) * z(i) + sig_1(i) * dz(i) + sig_2(i) * robot_state.dtheta(i);
    }
    // std::cout << "friction = " << friction.transpose() << std::endl;
    return friction;
}

Eigen::VectorXd Controller::friction_observer_PD(const RobotState & robot_state, Eigen::VectorXd& Control_input)  
{
    // error calculate
    Eigen::VectorXd de_nr = dtheta_n - robot_state.dtheta;
    Eigen::VectorXd e_nr = theta_n - robot_state.theta;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nq, nq);
    Eigen::VectorXd tau_f = -rotor_inertia_matrix_ * L * (de_nr + Lp*e_nr);
    Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( Control_input - robot_state.tau_J);   
    dtheta_n += step_time_ * ddtheta_n;
    theta_n += step_time_ * dtheta_n;
    tau_f_prev = tau_f;

    return tau_f;
}

Eigen::VectorXd Controller::friction_observer_L1_PD(const RobotState & robot_state, Eigen::VectorXd& Control_input)  
{
    // error calculate
    Eigen::VectorXd de_nr = dtheta_n - robot_state.dtheta;
    Eigen::VectorXd e_nr = theta_n - robot_state.theta;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nq, nq);
    Eigen::VectorXd sigma = -rotor_inertia_matrix_ * Gamma_ * (I + Gamma_ * step_time_).inverse()*((I + Gamma_p_ * step_time_)*de_nr + Gamma_p_*e_nr);
    // std::cout << "sigma = " << sigma.transpose() << std::endl;

    Eigen::VectorXd tau_f(nq);
    for (int i = 0; i < nq; ++i) {
       double alpha_i = std::exp(-K_lpf_(i) * step_time_);
       tau_f(i) = alpha_i * tau_f_prev(i) + (1.0 - alpha_i) * sigma(i);
    }
    // std::cout << "tau_f = " << tau_f.transpose() << std::endl;
    // state update
    Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( Control_input - robot_state.tau_J + sigma - tau_f );   
    dtheta_n += step_time_ * ddtheta_n;
    theta_n += step_time_ * dtheta_n;
    tau_f_prev = tau_f;

    return tau_f;
}

Eigen::VectorXd Controller::PD_controller_AS_GC(const RobotState & robot_state)
{
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity(); Kp.diagonal() << 50,50,50,50,50,50,50;
    Kd.setIdentity(); Kd.diagonal() << 3,3,3,3,3,3,3;

    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
    theta_des = robot_state.q_d + joint_stiffness_matrix_.inverse() * gravity_compensation_torque;

    Eigen::VectorXd control_motor_torque = -Kp * (robot_state.theta - theta_des)- Kd * (robot_state.dtheta - robot_state.dq_d)+ gravity_compensation_torque;
    return control_motor_torque;
}

// Eigen::VectorXd PD_controller_AS_GC_task_space(const RobotState & robot_state)
// {
//     Eigen::MatrixXd Kx(6,6),Dx(6,6),Kd(7,7);
//     Kx.setIdentity(); Kx.diagonal() << 1000,1000,1000,8,8,8;
//     Dx.setIdentity(); Dx.diagonal() << 100,100,100,2,2,2;
//     Kd.setIdentity(); Kd.diagonal() << 10,10,10,10,10,10,10;
    
//     Eigen::VectorXd x_err = robot_state.x - robot_state.x_d;    
//     Eigen::VectorXd dx_err = robot-state.dx - robot_state.dx_d;  
    
// }