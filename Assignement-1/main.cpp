#include <iostream>

#include <Eigen/Dense>

double deg_to_rad(double degrees)
{
    return degrees * 0.0174532925;
}

double rad_to_deg(double radians)
{
    return radians * 57.2957795;
}

// Equation (3.30) page 74, MR pre-print 2019
Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& vec) {
    Eigen::Matrix3d skew;
    skew << 0, -vec.z(), vec.y(),
            vec.z(), 0, -vec.x(),
            -vec.y(), vec.x(), 0;
    return skew;
}

// Equations 3.19 and 3.21 page 67, MR pre-print 2019
Eigen::Matrix3d rotation_matrix_from_frame_axes(const Eigen::Vector3d &x,
                                                const Eigen::Vector3d &y,
                                                const Eigen::Vector3d &z) {
    // Enforce right-handed coordinate system ergo the cross product of x and y must equal z
    Eigen::Vector3d xh = x.normalized();

    // Use Gram–Schmidt process to ensure yh is perpendicular to xh. https://en.wikipedia.org/wiki/Gram%E2%80%93Schmidt_process
    // This enforces R^T R = I
    Eigen::Vector3d yh = (y - xh * (xh.dot(y))).normalized();

    // Ensure right handedness with zh = xh X yh
    Eigen::Vector3d zh = xh.cross(yh).normalized();

    Eigen::Matrix3d rot_matrix;
        rot_matrix.col(0) = xh;
        rot_matrix.col(1) = yh;
        rot_matrix.col(2) = zh;
    return rot_matrix;
}

// Unnamed equation, page 72, MR pre-print 2019
Eigen::Matrix3d rotate_x(double degrees) {
    const double rad_angle = deg_to_rad(degrees);
    const double c = std::cos(rad_angle);
    const double s = std::sin(rad_angle);
    Eigen::Matrix3d rot_matrix;
        rot_matrix << 1, 0, 0,
                      0, c,-s,
                      0, s, c;
    return rot_matrix;
}

// Unnamed equation, page 72, MR pre-print 2019
Eigen::Matrix3d rotate_y(double degrees) {
    const double rad_angle = deg_to_rad(degrees);
    const double c = std::cos(rad_angle);
    const double s = std::sin(rad_angle);
    Eigen::Matrix3d rot_matrix;
    rot_matrix << c, 0, s,
                  0, 1, 0,
                 -s, 0, c;
    return rot_matrix;
}

// Unnamed equation, page 72, MR pre-print 2019
Eigen::Matrix3d rotate_z(double degrees) {
    const double rad_angle = deg_to_rad(degrees);
    const double c = std::cos(rad_angle);
    const double s = std::sin(rad_angle);
    Eigen::Matrix3d rot_matrix;
    rot_matrix << c,-s, 0,
                  s, c, 0,
                  0, 0, 1;
    return rot_matrix;
}

// Equation (3.51) page 101, MR pre-print 2019
// Rodrigues formula
Eigen::Matrix3d rotation_matrix_from_axis_angle(const Eigen::Vector3d &axis, double degrees){
    const double theta = deg_to_rad(degrees);

    const double c = std::cos(theta);
    const double s = std::sin(theta);

    const double n = axis.norm();

    const Eigen::Vector3d w = axis / n;
    Eigen::Matrix3d W = skew_symmetric(w);

    return Eigen::Matrix3d::Identity() + s * W + (1-c) * (W*W);
}

Eigen::Matrix3d rotation_matrix_from_euler_zyx(const Eigen::Vector3d &e){
        const double rot_x = e[2];
        const double rot_y = e[1];
        const double rot_z = e[0];

        Eigen::Matrix3d rot_matrix;
        rot_matrix = rotate_z(rot_z) * rotate_y(rot_y) * rotate_x(rot_x);
        return rot_matrix;
    }

// Equation (3.62) page 87, MR pre-print 2019
Eigen::Matrix4d transformation_matrix(const Eigen::Matrix3d &r, const Eigen::Vector3d &p){
    Eigen::Matrix4d transformation_matrix = Eigen::Matrix4d::Identity();
    transformation_matrix.topLeftCorner<3, 3>() = r;
    transformation_matrix.topRightCorner<3, 1>() = p;
    return transformation_matrix;
}

// Do not need to do translation since it is not in point space
void transform_vector() {
    const Eigen::Vector3d v_a (2.5, 3.0, -10.0);
    const Eigen::Vector3d Eigen_Rotations(60.0, 45.0, 0.0);

    const Eigen::Matrix3d rotation_matrix = rotation_matrix_from_euler_zyx(Eigen_Rotations);

    const Eigen::Vector3d v_w = rotation_matrix * v_a;

    std::cout << "v_w = [" << v_w.x() << ", " << v_w.y() << ", " << v_w.z() << "]\n";
}

/* Test functions below */
void skew_symmetric_test()
{
    Eigen::Matrix3d skew_matrix = skew_symmetric(Eigen::Vector3d{0.5, 0.5, 0.707107});
    std::cout << "Skew-symmetric matrix: " << std::endl;
    std::cout << skew_matrix << std::endl;
    std::cout << "Skew-symmetric matrix transposition: " << std::endl;
    std::cout << -skew_matrix.transpose() << std::endl;
}

void rotation_matrix_test()
{
    Eigen::Matrix3d rot = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45.0, -45.0, 90.0});
    Eigen::Matrix3d rot_aa = rotation_matrix_from_axis_angle(Eigen::Vector3d{0.8164966, 0.0,
    0.5773503}, 120.0);
    Eigen::Matrix3d rot_fa = rotation_matrix_from_frame_axes(Eigen::Vector3d{0.5, 0.5, 0.707107},
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
    Eigen::Matrix3d r = rotation_matrix_from_euler_zyx(Eigen::Vector3d{45, -45.0, 90.0});
    Eigen::Vector3d v{1.0, -2.0, 3.0};
    std::cout << "transformation_matrix: " << std::endl;
    std::cout << transformation_matrix(r, v) << std::endl;
}

int main()
{
    skew_symmetric_test();
    rotation_matrix_test();
    transformation_matrix_test();
    transform_vector();
    return 0;
}
