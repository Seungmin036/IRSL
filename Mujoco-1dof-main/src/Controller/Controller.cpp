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
    Gamma_ = Eigen::MatrixXd::Identity(nq, nq) * 10;
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
    {   Eigen::VectorXd lugre = lugre_friction(robot_state);
        Eigen::VectorXd tau_f_obs = friction_observer_PD(robot_state, u_prev);
        Eigen::VectorXd u_ctrl = this->PD_controller_AS_GC(robot_state) - lugre;// - tau_f_obs;
        u = u_ctrl;
        u_prev = u;
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
    Eigen::VectorXd dz(nq);

    for (int i = 0; i < nq; ++i)
    {
        double abs_dtheta = std::abs(robot_state.dtheta(i));
        double g = Fc(i) + (Fs(i) - Fc(i)) * std::exp(-std::pow(robot_state.dtheta(i) / vs(i), 2));

        if (abs_dtheta > 1e-6) 
            dz(i) = robot_state.dtheta(i) - (abs_dtheta / g) * z(i);
        else
            dz(i) = 0.0;

        z(i) += dz(i) * step_time_;

        friction(i) = sig_0(i) * z(i) + sig_1(i) * dz(i) + sig_2(i) * robot_state.dtheta(i);
    }
    std::cout << "friction = " << friction.transpose() << std::endl;
    return friction;
}

Eigen::VectorXd Controller::friction_observer_PD(const RobotState & robot_state, Eigen::VectorXd& Control_input)  
{
    Eigen::VectorXd ddtheta_n = rotor_inertia_matrix_.inverse() * ( Control_input - tau_j_prev + sigma_prev - tau_f_prev );   
    Eigen::VectorXd dtheta_n = dtheta_n_prev + step_time_ * ddtheta_n;
    for (int i = 0; i < nq; ++i) {
        if (std::abs(dtheta_n(i)) < 1e-4) {
            dtheta_n(i) = 0.0;
        }
    }
    Eigen::VectorXd theta_n = theta_n_prev + step_time_ * dtheta_n_prev;
    // error calculate
    Eigen::VectorXd de_nr = dtheta_n - robot_state.dtheta;
    Eigen::VectorXd e_nr = theta_n - robot_state.theta;
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(nq, nq);
    Eigen::VectorXd dde_nr = -(I + Gamma_ * step_time_).inverse() * Gamma_ * ((I + Gamma_p_ * step_time_) * de_nr + Gamma_p_ * e_nr);

    Eigen::VectorXd sigma = -rotor_inertia_matrix_ * Gamma_ * (step_time_ * dde_nr + (I + Gamma_p_ * step_time_) * de_nr + Gamma_p_ * e_nr);
    std::cout << "sigma = " << sigma.transpose() << std::endl;

    Eigen::VectorXd tau_f(nq);
    for (int i = 0; i < nq; ++i) {
       double alpha_i = std::exp(-K_lpf_(i) * step_time_);
       tau_f(i) = alpha_i * tau_f_prev(i) + (1.0 - alpha_i) * sigma(i);
    }
    // state update
    theta_n_prev = theta_n;
    dtheta_n_prev = dtheta_n;
    sigma_prev = sigma;
    tau_f_prev = tau_f;
    tau_j_prev = robot_state.tau_J;
    // 함수 내부에 static 파일 핸들 하나만 유지
    static std::ofstream file_friction_all("friction_full_log.csv");

    // 출력 포맷 지정
    file_friction_all << std::fixed << std::setprecision(6);

    // 1. sigma와 tau_f 교차 출력
    for (int i = 0; i < nq; ++i) {
        file_friction_all << sigma(i) << "," << tau_f(i) << ",";
    }

    // 2. error
    for (int i = 0; i < nq; ++i) {
        file_friction_all << e_nr(i) << ",";
    }

    // 3. de_nr
    for (int i = 0; i < nq; ++i) {
        file_friction_all << de_nr(i) << ",";
    }

    // 4. tau_j
    for (int i = 0; i < nq; ++i) {
        file_friction_all << robot_state.tau_J(i);
        if (i < nq - 1) file_friction_all << ",";  // 마지막 항목 이후에는 쉼표 안 붙임
    }

    // 줄바꿈
    file_friction_all << "\n";

    return sigma;
}

Eigen::VectorXd Controller::PD_controller_AS_GC(const RobotState & robot_state)
{
    Eigen::MatrixXd Kp(nq,nq),Kd(nv,nv);
    Kp.setIdentity();
    Kd.setIdentity();

    Kp.diagonal() << 50,50,50,50,50,50,50;
    Kd.diagonal() << 3,3,3,3,3,3,3;
    

    Eigen::VectorXd gravity_compensation_torque(nv), theta_des(nq);
    gravity_compensation_torque = robot.GetGravity(robot_state.q_d);
    theta_des = robot_state.q_d + joint_stiffness_matrix_.inverse() * gravity_compensation_torque;

    Eigen::VectorXd tau_c = -Kp * (robot_state.theta - theta_des)- Kd * (robot_state.dtheta - robot_state.dq_d)+ gravity_compensation_torque;
    Eigen::VectorXd control_motor_torque = tau_c;

    return control_motor_torque;
}