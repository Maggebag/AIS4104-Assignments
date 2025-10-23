#include <iostream>

#include "Assignment-1.hpp"

void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = robotics1::skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

void rotation_matrix_test()
{
    Eigen::Matrix3d rot = robotics1::rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa = robotics1::rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0,
                                                                            0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa = robotics1::rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
                                                                        Eigen::Vector3d{-0.5, -0.5, 0.707107}, Eigen::Vector3d{0.707107, -0.707107, 0.0});
    std::cout << "Rotation matrix from Euler: " << std::endl;
    std::cout << rot << std::endl << std::endl;
    std::cout << "Rotation matrix from axis-angle pair: " << std::endl;
    std::cout << rot_aa << std::endl << std::endl;
    std::cout << "Rotation matrix from frame axes: " << std::endl;
    std::cout << rot_fa << std::endl << std::endl;
}

void transformation_matrix_test()
{
    Eigen::Matrix3d r = robotics1::rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << robotics1::transformation_matrix(r, v) << std::endl;
}

int main()
{
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    robotics1::transform_vector();
    return 0;
}
