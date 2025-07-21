#include <Controller/ControllerNode.hpp>

ControllerNode::ControllerNode(ros::NodeHandle &nh, const pinocchio::Model &pinocchio_model) :
controller(pinocchio_model), nq(pinocchio_model.nq), nv(pinocchio_model.nv), nu(pinocchio_model.nv), is_rigid(false)
{  
    n=nh;
    joint_state_msg_pub = n.advertise<sensor_msgs::JointState>("joint_states",1000);
    
}
ControllerNode::~ControllerNode()
{

}

void ControllerNode::InitControllerNode(const mjModel* m, mjData* d)
{   
    robot_state_init.theta.resize(nq); 
    robot_state_init.dtheta.resize(nv);
    robot_state_init.q.resize(nq);
    robot_state_init.dq.resize(nv);

    robot_state_init.is_rigid = is_rigid;
    robot_state.is_rigid = is_rigid;

    Eigen::VectorXd theta_init(nq); theta_init << 0, 0, 0, 0, 0, 0, 0 ;
    // theta_init.setZero();

    if(is_rigid)
    {
        for(int i=0; i<nq; i++)
        {
            d->qpos[i] = theta_init(i);
        }
    }
    else
    {
        for(int i=0; i<nq; i++)
        {
            d->qpos[2*i] = theta_init(i);
        }
    }
    
    robot_state_init.theta = theta_init;
    robot_state_init.dtheta.setZero();
    robot_state_init.q = theta_init;
    robot_state_init.dq.setZero();

    //set desired value
    robot_state.q_d.resize(nq); robot_state.q_d << 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2;
    robot_state.dq_d.resize(nv); robot_state.dq_d.setZero();

    controller.InitController(robot_state_init);

    mj_forward(m,d);
}

void ControllerNode::UpdateRobotState(const mjModel* m, mjData* d, RobotState & robot_state)
{
    robot_state.theta.resize(nq); 
    robot_state.dtheta.resize(nv);
    robot_state.q.resize(nq);
    robot_state.dq.resize(nv);
    robot_state.tau_J.resize(nv);

    for(int i=0; i<nq; i++)
    {
        robot_state.theta(i) = d->qpos[2*i]; // motor side
        robot_state.dtheta(i)= d->qvel[2*i];
        robot_state.q(i)     = d->qpos[2*i] + d->qpos[2*i+1];
        robot_state.dq(i)    = d->qvel[2*i];
        robot_state.tau_J(i) = d->sensordata[3*i+2];
    }
}


void ControllerNode::Control_Loop(const mjModel* m, mjData* d)
{   
    //fill the data into robot_state structure
    UpdateRobotState(m,d,robot_state);
     
    Eigen::VectorXd u(nv);

    u = controller.GetControlInput(robot_state, CONTROLLER_SELECTOR::JOINT_PD_FRIC);

    //apply control input u
    for (int i=0; i<nv;i++)
    {
        d->ctrl[i] = u(i);
    }
    // std::cout << "Control_input = " << u.transpose() << std::endl;   
}
