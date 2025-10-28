#include <Eigen/Dense>
#include <functional>

#include "Assignment-3.hpp"
#include "Assignment-2.hpp"
#include "Assignment-1.hpp"

namespace robotics3 {

    Eigen::VectorXd std_vector_to_eigen(const std::vector<double>& v)
    {
        // This is not the most efficient way of doing this as we are just copying the information instead of referencing its memory
        // But who cares, sometimes copying it is more safe :)
        Eigen::VectorXd v_eigen(v.size());
        for (size_t i = 0; i < v.size(); ++i)
        {
            v_eigen(i) = v[i];
        }
        return v_eigen;
    }

    bool is_average_below_eps(const std::vector<double> &values, double eps, uint8_t n_values)
    {
        if (values.size() != n_values){return false;} //Not enough elements return false

        double sum = 0.0;
        for (size_t i = values.size() - n_values; i < values.size(); ++i)
        {
            sum += values[i];
        }

        const double avg = sum / n_values;
        return avg > eps;
    }

    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_space_chain()
    {
        // UR3e geometry constants
        double L1 = 0.24355;
        double L2 = 0.2132;
        double W1 = 0.13105;
        double W2 = 0.0921;
        double H1 = 0.1518;
        double H2 = 0.08535;

        // Axis directions
        std::vector<Eigen::Vector3d> w = {
            {0,0,1},
            {0,1,0},
            {0,1,0},
            {0,1,0},
            {0,0,-1},
            {0,1,0}
        };

        // Points on axes, know that screw axis are usually represented with cross product of the axis direction and a point on the axis
        // Which is what the v_i representation was in the ur3e_fk_screw function uses to find M.
        std::vector<Eigen::Vector3d> q = {
            {0.0,       0.0,        0.0      },
            {0.0,       0.0,        H1       },
            {0.0,       0.0,        H1       },
            {L1,        0.0,        H1       },
            {L1+L2,     0.0,        H1       },
            {L1+L2,     0.0,        H1-H2    }
        };

        std::vector<Eigen::VectorXd> s(6); // Create a vector of Eigen-Vectors to store our screw-axes

        for (int i = 0; i < s.size(); ++i) {
            s[i] = robotics2::screw_axis(q[i], w[i], 0.0);
        }

        Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
        M << -1, 0, 0, L1 + L2,
              0, 0, 1, W1 + W2,
              0, 1, 0, H1 - H2,
              0, 0, 0, 1;

        return {M, s};
    }

    Eigen::Matrix4d ur3e_space_fk(const Eigen::VectorXd &joint_positions)
    {
        auto [M,S]  = ur3e_space_chain();

        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

        for (int i = 0; i < S.size(); ++i) {
            const Eigen::Vector3d w = S[i].head<3>();
            const Eigen::Vector3d v = S[i].tail<3>();
            const double theta_deg = joint_positions[i] * robotics1::rad_to_deg;

            T *= robotics2::matrix_exponential(w, v, theta_deg);
        }

        return T*M;
    }

    // Create screw axis but from the frame of the end effector instead, "move backwards from the front"
    // Equation 4.16, page 146, MR pre-print 2019
    std::pair<Eigen::Matrix4d, std::vector<Eigen::VectorXd>> ur3e_body_chain()
    {
        auto [M,S]  = ur3e_space_chain(); // Get the space chain, we can use the adjoint to change into body frame
        const Eigen::Matrix4d Minv = M.inverse();
        const Eigen::Matrix<double, 6,6> AdMinv = robotics2::adjoint_matrix(Minv);

        std::vector<Eigen::VectorXd> B_screws(S.size(), Eigen::VectorXd::Zero(6));
        for (int i = 0; i < S.size(); ++i) {
            B_screws[i] = AdMinv * S[i];
        }

        return {M, B_screws};
    }

    Eigen::Matrix4d ur3e_body_fk(const Eigen::VectorXd &joint_positions)
    {
        auto [M,S]  = ur3e_body_chain();

        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

        for (int i = 0; i < S.size(); ++i) {
            const Eigen::Vector3d w = S[i].head<3>();
            const Eigen::Vector3d v = S[i].tail<3>();
            const double theta_deg = joint_positions[i] * robotics1::rad_to_deg;

            T *= robotics2::matrix_exponential(w, v, theta_deg);
        }

        return M*T;
    }

    // Section 6.2.1, page 225, MR pre-print 2019
    std::pair<uint32_t, double> newton_raphson_root_find(const std::function<double(double)> &f, double x_0, double dx_0, double eps)
    {
        const uint32_t max_iter = 1000;
        uint32_t iter = 0;
        double x = x_0;

        for (iter = 0; iter < max_iter; ++iter) {
            double fx = f(x);
            if (std::abs(fx) < eps) {
                // This would mean we have converged on a solution
                return {iter, x};
            }

            double dfdx = (f(x + dx_0) - fx) / dx_0;

            // Be defensive for division by zero
            if (std::abs(dfdx) < 1e-14) {
                break;
            }

            x -= fx / dfdx;
        }

        return {iter, x}; // Did not find a better convergence within iteration max
    }

    // Applied Gradient Decent in mathematical methods 3 so implementation is based on that python script
    std::pair<uint32_t, double> gradient_descent_minimize(const std::function<double(double)> &f, double x_0, double gamma, double dx_0, double eps)
    {
        constexpr uint32_t max_iter = 1000;
        uint32_t iter = 0;
        double x = x_0;

        auto derivative = [&](const double x_val) {
            return (f(x_val + dx_0) - f(x_val)) / dx_0;
        };

        for (iter = 0; iter < max_iter; ++iter) {
            double grad = derivative(x);
            if (std::abs(grad) < eps) {
                return {iter, grad}; //Convergence
            }

            x = x-gamma*grad;
        }

        return {iter, x}; // Max iterations
    }

    // Equation 5.11, page 178, MR pre-print 2019
    Eigen::MatrixXd ur3e_space_jacobian(const Eigen::VectorXd &current_joint_positions)
    {
        auto [_, S] = ur3e_space_chain();
        const int dof = S.size();

        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, dof);
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

        for (int i = 0; i < dof; ++i) {
            if (i == 0) {
                J.col(0) = S[0];
            } else {
                J.col(i) = robotics2::adjoint_matrix(T) * S[i];
            }

            const double theta_deg = current_joint_positions[i] * robotics1::rad_to_deg;
            const Eigen::Vector3d w = S[i].head<3>();
            const Eigen::Vector3d v = S[i].tail<3>();

            T = T * robotics2::matrix_exponential(w, v, theta_deg);
        }

        return J;
    }

    // Equation 5.18, page 183, MR pre-print 2019
    Eigen::MatrixXd ur3e_body_jacobian(const Eigen::VectorXd &current_joint_positions)
    {
        auto [_, B] = ur3e_body_chain();
        const int dof = B.size();

        Eigen::MatrixXd J = Eigen::MatrixXd::Zero(6, dof);
        Eigen::Matrix4d T = Eigen::Matrix4d::Identity();

        J.col(dof-1) = B[dof-1];

        for (int i = dof-2; i >= 0; --i)
        {
            const Eigen::Vector3d w = B[i+1].head<3>();
            const Eigen::Vector3d v = B[i+1].tail<3>();

            T = T * robotics2::matrix_exponential(w, v, -current_joint_positions[i+1] * robotics1::rad_to_deg);

            J.col(i) = robotics2::adjoint_matrix(T) * B[i];
        }

        return J;
    }

    std::pair<size_t, Eigen::VectorXd> ur3e_ik_body(const Eigen::Matrix4d &t_sd, const Eigen::VectorXd
    &current_joint_positions, double gamma, double v_e, double w_e)
    {

    }

}
