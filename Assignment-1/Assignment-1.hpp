
#ifndef AIS4104_ASSIGNMENT_N_ASSIGNMENT_1_HPP
#define AIS4104_ASSIGNMENT_N_ASSIGNMENT_1_HPP

#endif //AIS4104_ASSIGNMENT_N_ASSIGNMENT_1_HPP

#include <Eigen/Dense>

namespace robotics1 {

    double deg_to_rad(double degrees);
    double rad_to_deg(double radians);

    Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& vec);

    Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                    const Eigen::Vector3d &y,
                                                    const Eigen::Vector3d &z);

    Eigen::Matrix3d rotate_x(double degrees);
    Eigen::Matrix3d rotate_y(double degrees);
    Eigen::Matrix3d rotate_z(double degrees);

    Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis,
                                                    double degrees);

    Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e);

    Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r,
                                          const Eigen::Vector3d &p);

    void transform_vector();

} // namespace
