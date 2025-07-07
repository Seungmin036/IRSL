#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <vector>
#include <pinocchio/parsers/urdf.hpp>
#include <pinocchio/multibody/model.hpp>
#include <pinocchio/multibody/fwd.hpp>
#include <pinocchio/algorithm/rnea.hpp>
#include <pinocchio/algorithm/kinematics.hpp>
#include <pinocchio/algorithm/crba.hpp>
#include <pinocchio/algorithm/jacobian.hpp>
#include <pinocchio/algorithm/frames.hpp>

#include <data/CommonOperation.hpp>

namespace Dynamics
{
    class articulated_system
    {   private :
        Eigen::VectorXd gc_, gv_;
        pinocchio::Model model_;
        pinocchio::Data data_;

        public :
        
        articulated_system(const pinocchio::Model &model);

        void UpdateState(const Eigen::VectorXd& gc, const Eigen::VectorXd& gv);

        Eigen::MatrixXd GetMass(const Eigen::VectorXd &q);
        Eigen::MatrixXd GetCoriolis(const Eigen::VectorXd &q, const Eigen::VectorXd &dq);
        Eigen::MatrixXd GetGravity(const Eigen::VectorXd &q);
        Eigen::MatrixXd GetEEJacobianW(const Eigen::VectorXd &q, const std::string & ee_name);
        Eigen::Affine3d GetEEpose(const Eigen::VectorXd &q, const std::string & ee_name);

        Eigen::VectorXd GetBaseWrench_RNEA(const Eigen::VectorXd &gacc);

    };
}