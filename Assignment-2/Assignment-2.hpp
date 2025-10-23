
#ifndef AIS4104_ASSIGNMENT_N_ASSIGNMENT_2_HPP
#define AIS4104_ASSIGNMENT_N_ASSIGNMENT_2_HPP

#endif //AIS4104_ASSIGNMENT_N_ASSIGNMENT_2_HPP

#include <Eigen/Dense>

namespace robotics2 {

    Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r);

    Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v);

    Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h);

    Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf);

    double cot(double x);

    Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d& e_deg);

    Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta);

    std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r);

    Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta);

    std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t);

    void print_pose(const std::string &label, const Eigen::Matrix4d &tf);

    void task2_a();

    void task2_b();

    Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions);

    Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions);

    void test_3r_planar_tf();

    void test_3r_planar_screw();

    Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions);

    Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions);

    void test_ur3_fk_screw();

    void test_ur3_fk_transform();


}