#include <iostream>

#include "Assignment-3.hpp"
#include "Assignment-2.hpp"
#include "Assignment-1.hpp"

void ur3e_test_fk()
{
    std::cout << "Forward kinematics tests" << std::endl;
    robotics2::print_pose("Space-1", robotics3::ur3e_space_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    robotics2::print_pose("Body-1", robotics3::ur3e_body_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    std::cout << std::endl;
    robotics2::print_pose("Space-2", robotics3::ur3e_space_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    robotics2::print_pose("Body-2", robotics3::ur3e_body_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, -90.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    std::cout << std::endl;
    robotics2::print_pose("Space-3", robotics3::ur3e_space_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    robotics2::print_pose("Body-3" ,robotics3::ur3e_body_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -180.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    std::cout << std::endl;
    robotics2::print_pose("Space-4", robotics3::ur3e_space_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
    robotics2::print_pose("Body-4", robotics3::ur3e_body_fk(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -90.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad));
};

void test_newton_raphson_root_find(const std::function<double(double)> &f, double x0)
{
    auto [iterations, x_hat] = robotics3::newton_raphson_root_find(f, x0);
    std::cout << "NR root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" <<
    f(x_hat) << std::endl;
}

void test_gradient_descent_minimize(const std::function<double(double)> &f, double x0)
{
    auto [iterations, x_hat] = robotics3::gradient_descent_minimize(f, x0);
    std::cout << "GD root f, x0=" << x0 << " -> it=" << iterations << " x=" << x_hat << " f(x)=" <<
    f(x_hat) << std::endl;
}

void test_optimizations()
{
    std::cout << "Root finding tests" << std::endl;
    auto f1 = [](double x)
    {
        return (x - 3.0) * (x - 3.0) - 1.0;
    };
    test_newton_raphson_root_find(f1, -20.0);
    test_gradient_descent_minimize(f1, -20.0);

    // Run another test with a different starting point
    test_newton_raphson_root_find(f1, 20.0);
    test_gradient_descent_minimize(f1, 20.0);
}

void ur3e_test_jacobian(const Eigen::VectorXd &joint_positions)
{
    Eigen::Matrix4d tsb = robotics3::ur3e_body_fk(joint_positions);
    auto [m, space_screws] = robotics3::ur3e_space_chain();
    Eigen::MatrixXd jb = robotics3::ur3e_body_jacobian(joint_positions);
    Eigen::MatrixXd js = robotics3::ur3e_space_jacobian(joint_positions);
    Eigen::MatrixXd ad_tsb = robotics2::adjoint_matrix(tsb);
    Eigen::MatrixXd ad_tbs = robotics2::adjoint_matrix(tsb.inverse());
    std::cout << "Jb: " << std::endl << jb << std::endl << "Ad_tbs*Js:" << std::endl << ad_tbs * js <<
    std::endl << std::endl;
    std::cout << "Js: " << std::endl << js << std::endl << "Ad_tsb*Jb:" << std::endl << ad_tsb * jb <<
    std::endl << std::endl;
    std::cout << "d Jb: " << std::endl << jb - ad_tbs * js << std::endl << std::endl;
    std::cout << "d Js: " << std::endl << js - ad_tsb * jb << std::endl << std::endl;
}

void ur3e_test_jacobian()
{
    std::cout << "Jacobian matrix tests" << std::endl;
    ur3e_test_jacobian(robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad);
    ur3e_test_jacobian(robotics3::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0, -50.0}) *
    robotics1::deg_to_rad);
}

void ur3e_ik_test_pose(const Eigen::Vector3d &pos, const Eigen::Vector3d &zyx, const Eigen::VectorXd
&j0)
{
    std::cout << "Test from pose" << std::endl;
    Eigen::Matrix4d t_sd = robotics1::transformation_matrix(robotics1::rotation_matrix_from_euler_zyx(zyx * robotics1::rad_to_deg), pos);
    auto [iterations, j_ik] = robotics3::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = robotics3::ur3e_body_fk(j_ik);
    robotics2::print_pose(" IK pose", t_ik);
    robotics2::print_pose("Desired pose", t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * robotics1::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * robotics1::rad_to_deg << std::endl << std::endl;
}

void ur3e_ik_test_configuration(const Eigen::VectorXd &joint_positions, const Eigen::VectorXd &j0)
{
    std::cout << "Test from configuration" << std::endl;
    Eigen::Matrix4d t_sd = robotics3::ur3e_space_fk(joint_positions);
    auto [iterations, j_ik] = robotics3::ur3e_ik_body(t_sd, j0);
    Eigen::Matrix4d t_ik = robotics3::ur3e_body_fk(j_ik);
    robotics2::print_pose("IK pose", t_ik);
    robotics2::print_pose("Desired pose", t_sd);
    std::cout << "Converged after " << iterations << " iterations" << std::endl;
    std::cout << "J_0: " << j0.transpose() * robotics1::rad_to_deg << std::endl;
    std::cout << "J_d: " << joint_positions.transpose() * robotics1::rad_to_deg << std::endl;
    std::cout << "J_ik: " << j_ik.transpose() * robotics1::rad_to_deg << std::endl << std::endl;
}

void ur3e_ik_test()
{
    std::cout << "----------------------------------------" << std::endl;
    Eigen::VectorXd j_t0 = robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, 0.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad;
    Eigen::VectorXd j_t1 = robotics3::std_vector_to_eigen(std::vector<double>{0.0, 0.0, -89.0, 0.0, 0.0, 0.0}) *
    robotics1::deg_to_rad;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} *
    robotics1::deg_to_rad, j_t0);
    std::cout << "----------------------------------------" << std::endl;
    ur3e_ik_test_pose(Eigen::Vector3d{0.3289, 0.22315, 0.36505}, Eigen::Vector3d{0.0, 90.0, -90.0} *
    robotics1::deg_to_rad, j_t1);
    std::cout << "----------------------------------------" << std::endl;
    Eigen::VectorXd j_t2 = robotics3::std_vector_to_eigen(std::vector<double>{50.0, -30.0, 20, 0.0, -30.0, 50.0})
    * robotics1::deg_to_rad;
    Eigen::VectorXd j_d1 = robotics3::std_vector_to_eigen(std::vector<double>{45.0, -20.0, 10.0, 2.5, 30.0,
    -50.0}) * robotics1::deg_to_rad;
    ur3e_ik_test_configuration(j_d1, j_t0);
    std::cout << "----------------------------------------" << std::endl;
    ur3e_ik_test_configuration(j_d1, j_t2);
}

int main()
{
    ur3e_test_fk();
    test_optimizations();
    // We find the correct roots of f(x) ≈ 0 at x = 2 and x = 4 using Newton–Raphson.
    // Gradient descent does not find the roots, because it is designed to minimize the
    // function value (i.e., solve f'(x) = 0) rather than find where f(x) = 0.
    ur3e_test_jacobian();
    ur3e_ik_test();
    // IK solver works well but does not necesarily lock on to the solution we want. For example for the pose testing it finds another correct pose that aquires the same TCP, but not necesarily the correct orientation. This could be solved with a bit more advanced code and some cost functions where we penalise large changes in orientation
    // Also same goes for the test configuration. We find the correct pose and orientation but the jacobian is not exactly at the correct solution. We know that multiple joint configurations can realize the same TCP pose, and this is what we see an example of here.
    return 0;
}