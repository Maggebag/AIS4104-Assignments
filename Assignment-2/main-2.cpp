#include <iostream>
#include <numbers>

#include "Assignment-2.hpp"
#include "Assignment-1.hpp"

// My own test functions since none where given
void test_Euler_zyx()
{
    constexpr double yaw = std::numbers::pi / 4;
    constexpr double pitch = - std::numbers::pi / 2;
    constexpr double roll = std::numbers::pi / 8;

    // Why did I make it take degrees input, in hindsight makes it really annoying to work with
    const Eigen::Vector3d euler_zyx_deg(robotics1::rad_to_deg(yaw), robotics1::rad_to_deg(pitch), robotics1::rad_to_deg(roll));

    const Eigen::Matrix3d R = robotics1::rotation_matrix_from_euler_zyx(euler_zyx_deg);

    const Eigen::Vector3d euler_zyx(yaw, pitch, roll);

    std::cout << "Original ZYX: " << euler_zyx.transpose() << std::endl;

    Eigen::Vector3d euler_zyx_recovered = robotics2::euler_zyx_from_rotation_matrix(R);

    std::cout << "Recovered ZYX: " << euler_zyx_recovered.transpose() << "\n" << std::endl;
}

void test_twist_func()
{
    const Eigen::Vector3d w(10,20,30);
    const Eigen::Vector3d v(15,7,0);

    std::cout << "Original v: " << v.transpose() << " | Original w: " << w.transpose() << std::endl;
    std::cout << "New Twist Vector: " << robotics2::twist(w, v).transpose() << "\n" << std::endl;
}

void test_screw_func()
{
    Eigen::Vector3d s(0, 0, 1);

    Eigen::Vector3d q(1, 0, 0);

    double h = 0.5;

    Eigen::VectorXd S = robotics2::screw_axis(q, s, h);

    std::cout << "Screw axis S = " << S.transpose() << "\n" << std::endl;
}

void test_adjoint_matrix()
{
    Eigen::Matrix3d R = robotics1::rotate_z(45);

    Eigen::Vector3d p (1,2,3);

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3,3>() = R;
    T.topRightCorner<3,1>() = p;

    Eigen::MatrixXd Adj = robotics2::adjoint_matrix(T);

    std::cout << "AdjointMatrix: \n " << Adj << "\n" << std::endl;
}

// Could make test functions for further functions in task 3, however have not done it. I trust my code... famous last words

int main()
{
    /*test_Euler_zyx();
    test_twist_func();
    test_screw_func();
    test_adjoint_matrix(); */
    robotics2::task2_a();
    robotics2::task2_b();
    robotics2::test_3r_planar_tf();
    robotics2::test_3r_planar_screw();
    robotics2::test_ur3_fk_screw();
    robotics2::test_ur3_fk_transform();
    return 0;
}