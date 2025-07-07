#include "Dynamics/Dynamics.hpp"

namespace Dynamics
{   
    articulated_system::articulated_system(const pinocchio::Model &model) : data_(model)
    {   
        model_ = model;
    }

    void articulated_system::UpdateState(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv)
    {   
        gc_ = gc; 
        gv_ = gv; 
        pinocchio::forwardKinematics(model_, data_, gc, gv);
    }

    Eigen::VectorXd articulated_system::GetBaseWrench_RNEA(const Eigen::VectorXd &gacc)
    {
        pinocchio::rnea(model_, data_, gc_, gv_, gacc);
        std::vector<int> root_index;

        for(int i=0; i<model_.njoints; i++)
        {
            if(model_.parents[i]==0 && i!=0)
            {
                root_index.push_back(i);
            }
        }

        std::vector<Eigen::VectorXd> root_joint_wrenchs;

        for(int i=0; i<root_index.size(); i++)
        {   
            int root_joint_index = root_index[i];
            Eigen::VectorXd root_joint_wrench(6); root_joint_wrench.setZero();
            auto SE3 = data_.oMi[root_joint_index];
            auto f = data_.f[root_joint_index];
            root_joint_wrench.head(3) = SE3.rotation()*f.linear();
            root_joint_wrench.tail(3) = SE3.rotation()*f.angular();

            root_joint_wrenchs.push_back(root_joint_wrench);
        }

        Eigen::Vector3d base_poscom_w;
        auto base_SE3 = data_.oMi[0];
        base_poscom_w = base_SE3.translation() + base_SE3.rotation()*model_.inertias[0].lever();
        Eigen::VectorXd base_wrench(6); base_wrench.setZero();
        

        for(int i=0; i<root_index.size(); i++)
        {   
            Eigen::VectorXd root_joint_wrench = root_joint_wrenchs[i];
            int root_joint_index = root_index[i];
            
            auto root_SE3 = data_.oMi[root_joint_index];

            base_wrench.head(3) += root_joint_wrench.head(3);
            base_wrench.tail(3) += root_joint_wrench.tail(3) 
                                + skew(root_SE3.translation() - base_SE3.translation())*root_joint_wrench.head(3) ;
                                
        }
        
        base_wrench.head(3) += - model_.inertias[0].mass()*model_.gravity.linear();
        base_wrench.tail(3) += - skew(base_poscom_w - base_SE3.translation())*model_.inertias[0].mass()*model_.gravity.linear();

        return base_wrench;
    }

    Eigen::MatrixXd articulated_system::GetMass(const Eigen::VectorXd &q)
    {
        pinocchio::crba(model_,data_,q); // compute mass matrix
    
        Eigen::MatrixXd temp = data_.M;
        temp += temp.transpose().eval();
        temp.diagonal() = 0.5*temp.diagonal();

        return temp;
    }

    Eigen::MatrixXd articulated_system::GetCoriolis(const Eigen::VectorXd &q, const Eigen::VectorXd &dq)
    {
        pinocchio::computeCoriolisMatrix(model_,data_, q, dq); // compute Coriolis matrix
        return data_.C;
    }

    Eigen::MatrixXd articulated_system::GetGravity(const Eigen::VectorXd &q)
    {
        pinocchio::computeGeneralizedGravity(model_, data_, q);
        return data_.g;
    }

    Eigen::MatrixXd articulated_system::GetEEJacobianW(const Eigen::VectorXd &q, const std::string & ee_name)
    {   
        pinocchio::framesForwardKinematics(model_, data_, q); // compute FK
        // bool is_exist = model_.existFrame(end_effector_name,pinocchio::FrameType::FIXED_JOINT);
        auto ee_id = model_.getFrameId(ee_name,pinocchio::FrameType::FIXED_JOINT);
        auto SE3 = data_.oMf[ee_id]; // 1~7 joint frame index, 8 : end-effector index

        Eigen::MatrixXd Jacobian_b(6,model_.nv); Jacobian_b.setZero();
        pinocchio::frameJacobian(model_, data_, q, ee_id, Jacobian_b);

        Eigen::Matrix<double,6,6> gen_rot; gen_rot.setZero();
        Eigen::MatrixXd Jacobian_w(6,model_.nv); Jacobian_w.setZero();
        gen_rot.topLeftCorner(3,3) = SE3.rotation();
        gen_rot.bottomRightCorner(3,3) = SE3.rotation(); 

        Jacobian_w = gen_rot * Jacobian_b;

        return Jacobian_w;   
    }

    Eigen::Affine3d articulated_system::GetEEpose(const Eigen::VectorXd &q, const std::string & ee_name)
    {
        Eigen::Affine3d ee_htm;
        pinocchio::framesForwardKinematics(model_, data_, q); // compute FK
        // bool is_exist = model_.existFrame(end_effector_name,pinocchio::FrameType::FIXED_JOINT);
        auto ee_id = model_.getFrameId(ee_name,pinocchio::FrameType::FIXED_JOINT);
        auto SE3 = data_.oMf[ee_id]; // 1~7 joint frame index, 8 : end-effector index
        ee_htm = SE3.toHomogeneousMatrix();

        return ee_htm;
    }

}