
#pragma once

#include <Eigen/Dense>
#include <functional>
#include <numbers>

namespace robotics3 {

    Eigen::VectorXd std_vector_to_eigen(const std::vector<double> &v);

    bool is_average_below_eps(const std::vector<double> &values, double eps = 10e-7, uint8_t n_values = 5u);

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_space_chain();

    Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions);

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain();

    Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions);

    std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double
    x_0, double dx_0 = 0.5, double eps = 10e-7);

    std::pair<uint32_t, double> gradient_descent_minimize(const std::function<double(double)> &f, double
    x_0, double gamma = 0.1, double dx_0 = 0.5, double eps = 10e-7);

    Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions);

    Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions);

    std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd
    &current_joint_positions, double gamma = 1e-2, double v_e = 2e-3, double w_e = 2e-3);

}