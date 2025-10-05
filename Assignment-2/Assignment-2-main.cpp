#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>

#include <Eigen/Dense>

double deg_to_rad(double degrees)
{
    return degrees * 0.0174532925;
}

double rad_to_deg(double radians)
{
    return radians * 57.2957795;
}

// Utils from Assignment 1

// Equation (3.30) page 74, MR pre-print 2019
Eigen::Matrix3d skew_symmetric(const Eigen::Vector3d& vec) {
    Eigen::Matrix3d skew;
    skew << 0, -vec.z(), vec.y(),
            vec.z(), 0, -vec.x(),
            -vec.y(), vec.x(), 0;
    return skew;
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

// Euler rotation is just a way of creating a rotation matrix with a set order.
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

/* ---------------- Here comes the assignments ----------------*/
/* ----------------          Task 1           ---------------- */

// Equations from MR-preprint 2019 Appendix B chapter B.1.1
Eigen::Vector3d euler_zyx_from_rotation_matrix(const Eigen::Matrix3d &r)
{
    Eigen::Vector3d euler_zyx; // Init the vector for returning the rotations

    double alpha = atan2(r(1, 0), r(0, 0));
    double beta = atan2(-r(2,0), sqrt(std::pow(r(0,0),2) + std::pow(r(1,0),2)));
    double gamma = atan2(r(2,1), r(2,2));

    euler_zyx << alpha, beta, gamma;

    return euler_zyx;
}

// Equation 3.70 page 96, MR pre-print 2019
Eigen::VectorXd twist(const Eigen::Vector3d &w, const Eigen::Vector3d &v)
{
    Eigen::VectorXd twist(6);
    twist << w, v;
    return twist;
}

// Equation on page 101, MR pre-print 2019
Eigen::VectorXd screw_axis(const Eigen::Vector3d &q, const Eigen::Vector3d &s, double h)
{
    Eigen::Vector3d v = -s.cross(q) + h * s;
    Eigen::VectorXd axis(6);
    axis << s, v;
    return axis;
}

// Definition 3.20 on page 98, MR pre-print 2019
Eigen::MatrixXd adjoint_matrix(const Eigen::Matrix4d &tf)
{
    const Eigen::Matrix3d R = tf.topLeftCorner<3,3>();
    const Eigen::Vector3d p = tf.topRightCorner<3,1>();

    Eigen::Matrix<double,6,6> Adj;
    Adj.setZero();
    Adj.topLeftCorner<3,3>() = R;
    Adj.bottomLeftCorner<3,3>() = skew_symmetric(p) * R;
    Adj.bottomRightCorner<3,3>() = R;

    return Adj;
}

double cot(double x)
{
    return 1 / (std::sin(x) / std::cos(x)); // Use the known property that tan = sin/cos
}

/* ---------------- Task 2 ---------------- */
Eigen::Matrix3d rotation_matrix_from_euler_yzx(const Eigen::Vector3d& e_deg) {
    return rotate_y(e_deg[0]) * rotate_z(e_deg[1]) * rotate_x(e_deg[2]);
}

void task2_a()
{
    Eigen::Vector3d f_w (-30, 0, 0);
    Eigen::Vector3d m_s (0, 0, 2);
    Eigen::Vector3d e_ws(60, -60, 0);

    // Rotation from sensor to world
    Eigen::Matrix3d R_ws = rotation_matrix_from_euler_yzx(e_ws);

    Eigen::Vector3d m_w = R_ws * m_s;

    Eigen::Vector3d f_s = R_ws.transpose() * f_w;

    std::cout << "f_w: " << f_w.transpose() << "\n" << std::endl;
    std::cout << "m_w: " << m_w.transpose() << "\n" << std::endl;
    std::cout << "f_s: " << f_s.transpose() << "\n" << std::endl;
    std::cout << "m_s: " << m_s.transpose() << "\n" << std::endl;
}

void task2_b()
{
    // Define variables from example
    const double m_apple = 0.1; //kg
    const double g = 10; // m/s^2
    const double m_hand = 0.5; //kg
    const double L1 = 0.1; // m
    const double L2 = 0.15; // m

    Eigen::Vector3d w_h (0,0,0);
    Eigen::Vector3d f_h (0, -g*m_hand, 0);

    Eigen::Vector3d w_a (0, 0, 0);
    Eigen::Vector3d f_a (0, 0, g*m_apple);

    Eigen::VectorXd F_h = twist(w_h, f_h);
    Eigen::VectorXd F_a = twist(w_a, f_a);

    Eigen::Matrix3d R_h = Eigen::Matrix3d::Identity();
    Eigen::Matrix3d R_a;
    R_a << 1, 0, 0,
    0, 0, 1,
    0, -1, 0;

    Eigen::Vector3d p_h (-L1, 0, 0);
    Eigen::Vector3d p_a (-(L1+L2), 0, 0);

    Eigen::Matrix4d T_hf = transformation_matrix(R_h, p_h);
    Eigen::Matrix4d T_af = transformation_matrix(R_a, p_a);

    Eigen::MatrixXd Adj_hf = adjoint_matrix(T_hf);
    Eigen::MatrixXd Adj_af = adjoint_matrix(T_af);

    Eigen::VectorXd F_f = Adj_hf.transpose() * F_h + Adj_af.transpose()*F_a;

    std::cout << "F_f: " << F_f.transpose() << "\n" << std::endl;
}

/* ---------------- Task 3 ---------------- */
// Equation 3.51 page 82, MR pre-print 2019
Eigen::Matrix3d matrix_exponential(const Eigen::Vector3d &w, double theta)
{
    double rad = deg_to_rad(theta);
    Eigen::Matrix3d R = Eigen::Matrix3d::Identity() + sin(rad) * skew_symmetric(w) + (1-cos(rad)) * skew_symmetric(w) * skew_symmetric(w);
    return R;
}

// Equations 3.58 - 3.61, pages 85 - 86, MR pre-print 2019
std::pair<Eigen::Vector3d, double> matrix_logarithm(const Eigen::Matrix3d &r)
{
    double theta;
    Eigen::Vector3d w;

    if (r == Eigen::Matrix3d::Identity()) {
        theta = 0;
    }
    else {
        double trace_r = r.trace(); // Built-in eigen function for trace
        if (trace_r == -1) {
            theta = M_PI;
            w = (1 / sqrt(2 * (1 + r(2,2)))) * Eigen::Vector3d (r(0,2), r(1,2), 1 + r(2,2));
        }
        else {
            theta = acos(0.5 * (trace_r - 1));
            double w_n = 1 / 2*sin(theta);

            double w_1 = w_n * (r(2, 1) - r(1, 2));
            double w_2 = w_n * (r(0, 2) - r(2, 0));
            double w_3 = w_n * (r(1, 0) - r(0, 1));

            w = Eigen::Vector3d(w_1, w_2, w_3);
        }
    }

    return std::pair<Eigen::Vector3d, double>(w, theta);
}

// Preposition 3.25 on page 103, MR pre-print 2019
Eigen::Matrix4d matrix_exponential(const Eigen::Vector3d &w, const Eigen::Vector3d &v, double theta)
{
    double rad = deg_to_rad(theta);

    Eigen::Matrix3d R = matrix_exponential(w, theta);
    Eigen::Vector3d v_theta = ((Eigen::Matrix3d::Identity() * rad) + ((1-cos(rad)) * skew_symmetric(w)) + (rad - sin(rad)) * skew_symmetric(w) * skew_symmetric(w)) * v;

    Eigen::Matrix4d T = transformation_matrix(R, v_theta);
    return T;
}

// Algorithm on page 104, MR pre-print 2019
std::pair<Eigen::VectorXd, double> matrix_logarithm(const Eigen::Matrix4d &t)
{
    const Eigen::Matrix3d R = t.topLeftCorner<3,3>();
    const Eigen::Vector3d p = t.topRightCorner<3,1>();

    Eigen::Vector3d w;
    Eigen::Vector3d v;
    double theta;

    if (R == Eigen::Matrix3d::Identity()) {
        w = Eigen::Vector3d::Zero();
        v = p / p.norm();
        theta = p.norm();
    }
    else {
        std::pair<Eigen::Vector3d, double> m_log = matrix_logarithm(R);
        w = m_log.first;
        theta = m_log.second;
        const Eigen::Matrix3d skew_w = skew_symmetric(w);
        v = (((1/theta) * Eigen::Matrix3d::Identity()) - 0.5 * skew_w + ((1/theta) - 0.5* cot(theta/2)) * skew_w*skew_w) * p;
    }

    return std::pair<Eigen::VectorXd, double>(twist(w,v), theta);
}

/* ---------------- Task 4 ---------------- */
void print_pose(const std::string &label, const Eigen::Matrix4d &tf)
{
    Eigen::Matrix3d R = tf.topLeftCorner<3,3>();
    Eigen::Vector3d p = tf.topRightCorner<3,1>();
    Eigen::Vector3d e_ZYX = euler_zyx_from_rotation_matrix(R);
    std::cout << "Pose: \n" << "Euler ZYX: \n" << e_ZYX << "\n" << "Position: " << p.transpose() << "\n" << std::endl;
}

Eigen::Matrix4d planar_3r_fk_transform(const std::vector<double> &joint_positions)
{
    const double L1 = 10;
    const double L2 = 10;
    const double L3 = 10;

    Eigen::Matrix4d T_01 = transformation_matrix(rotate_z(joint_positions[0]), Eigen::Vector3d::Zero());
    Eigen::Matrix4d T_12 = transformation_matrix(rotate_z(joint_positions[1]), Eigen::Vector3d (L1,0,0));
    Eigen::Matrix4d T_23 = transformation_matrix(rotate_z(joint_positions[2]), Eigen::Vector3d (L2,0,0));
    Eigen::Matrix4d T_34 = transformation_matrix(Eigen::Matrix3d::Identity(), Eigen::Vector3d (L3,0,0));

    return T_01 * T_12 * T_23 * T_34;
}

// Example 4.2 on page 142, MR pre-print 2019
Eigen::Matrix4d planar_3r_fk_screw(const std::vector<double> &joint_positions)
{
    const double L1 = 10;
    const double L2 = 10;
    const double L3 = 10;

    Eigen::Vector3d w1 (0,0,1);
    Eigen::Vector3d v1 (0,0,0);
    Eigen::Vector3d w2 (0,0,1);
    Eigen::Vector3d v2 (0, -L1, 0);
    Eigen::Vector3d w3 (0,0, 1);
    Eigen::Vector3d v3 (0, -(L1 + L2), 0);

    Eigen::Matrix4d M = Eigen::Matrix4d::Identity();
    M(0,3) = L1+L2+L3;

    Eigen::Matrix4d e_1 = matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Matrix4d e_2 = matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Matrix4d e_3 = matrix_exponential(w3, v3, joint_positions[2]);

    return e_1 * e_2 * e_3 * M;
}

void test_3r_planar_tf()
{
    const std::vector<std::vector<double>> joint_positions = {
        {  0.0,   0.0,   0.0 },   // j1
        { 90.0,   0.0,   0.0 },   // j2
        {  0.0,  90.0,   0.0 },   // j3
        {  0.0,   0.0,  90.0 },   // j4
        { 10.0, -15.0,   2.75 }   // j5
    };
    std::cout << "Planar 3r transform: \n" << std::endl;
    for (int i = 0; i < joint_positions.size(); ++i ) {
        Eigen::Matrix4d T = planar_3r_fk_transform(joint_positions[i]);
        print_pose("T", T);
    }
}

void test_3r_planar_screw()
{
    const std::vector<std::vector<double>> joint_positions = {
        {  0.0,   0.0,   0.0 },   // j1
        { 90.0,   0.0,   0.0 },   // j2
        {  0.0,  90.0,   0.0 },   // j3
        {  0.0,   0.0,  90.0 },   // j4
        { 10.0, -15.0,   2.75 }   // j5
    };
    std::cout << "Planar 3r transform using POE: \n" << std::endl;
    for (int i = 0; i < joint_positions.size(); ++i ) {
        Eigen::Matrix4d T = planar_3r_fk_screw(joint_positions[i]);
        print_pose("T", T);
    }
}

/* ---------------- Task 5 ---------------- */
// Example 4.5 on page 145, MR pre-print 2019
// Variables changed for ur3e
Eigen::Matrix4d ur3e_fk_screw(const std::vector<double> &joint_positions)
{
    double L1 = 0.24355;
    double L2 = 0.2132;
    double W1 = 0.13105;
    double W2 = 0.0921;
    double H1 = 0.1518;
    double H2 = 0.08535;

    Eigen::Vector3d w1 (0,0,1);
    Eigen::Vector3d v1 (0, 0, 0);
    Eigen::Vector3d w2 (0, 1, 0);
    Eigen::Vector3d v2 (-H1, 0, 0);
    Eigen::Vector3d w3 (0, 1, 0);
    Eigen::Vector3d v3 (-H1, 0, L1);
    Eigen::Vector3d w4 (0, 1, 0);
    Eigen::Vector3d v4 (-H1, 0, L1 + L2);
    Eigen::Vector3d w5 (0, 0, -1);
    Eigen::Vector3d v5 (-W1, L1 + L2, 0);
    Eigen::Vector3d w6 (0, 1, 0);
    Eigen::Vector3d v6 (H2 - H1, 0, L1 + L2);

    Eigen::Matrix4d M { {-1, 0, 0, L1 + L2},
                            {0, 0, 1, W1 + W2 },
                            {0, 1, 0, H1 - H2 },
                            {0, 0, 0, 1 }};

    Eigen::Matrix4d e_1 = matrix_exponential(w1, v1, joint_positions[0]);
    Eigen::Matrix4d e_2 = matrix_exponential(w2, v2, joint_positions[1]);
    Eigen::Matrix4d e_3 = matrix_exponential(w3, v3, joint_positions[2]);
    Eigen::Matrix4d e_4 = matrix_exponential(w4, v4, joint_positions[3]);
    Eigen::Matrix4d e_5 = matrix_exponential(w5, v5, joint_positions[4]);
    Eigen::Matrix4d e_6 = matrix_exponential(w6, v6, joint_positions[5]);

    return e_1 * e_2 * e_3 * e_4  * e_5 * e_6 * M;
}

Eigen::Matrix4d ur3e_fk_transform(const std::vector<double> &joint_positions)
{
    double L1 = 0.24355;
    double L2 = 0.2132;
    double W1 = 0.13105;
    double W2 = 0.0921;
    double H1 = 0.1518;
    double H2 = 0.08535;

    const double q1 = deg_to_rad(joint_positions[0]);
    const double q2 = deg_to_rad(joint_positions[1]);
    const double q3 = deg_to_rad(joint_positions[2]);
    const double q4 = deg_to_rad(joint_positions[3]);
    const double q5 = deg_to_rad(joint_positions[4]);
    const double q6 = deg_to_rad(joint_positions[5]);

    // Can create transformation matrices directly from DH table
    // Had this function from an assignment from the bachelor studies
    auto DH = [](double th, double a, double d, double al) {
        const double ct = std::cos(th),  st = std::sin(th);
        const double ca = std::cos(al),  sa = std::sin(al);
        Eigen::Matrix4d T;
        T <<  ct, -st*ca,  st*sa,  a*ct,
              st,  ct*ca, -ct*sa,  a*st,
               0,      sa,     ca,     d,
               0,       0,      0,     1;
        return T;
    };

    // Create transformation matrices directly from DH-parameters (standard DH)
    Eigen::Matrix4d T_01 = DH(q1,      0.0,  H1,  M_PI/2);
    Eigen::Matrix4d T_12 = DH(q2,    -L1,   0.0,  0.0    );
    Eigen::Matrix4d T_23 = DH(q3,    -L2,   0.0,  0.0    );
    Eigen::Matrix4d T_34 = DH(q4,     0.0,  W1,  M_PI/2 );
    Eigen::Matrix4d T_45 = DH(q5,     0.0,  H2, -M_PI/2 );
    Eigen::Matrix4d T_56 = DH(q6,     0.0,  W2,  0.0    );

    return T_01 * T_12 * T_23 * T_34 * T_45 * T_56;
}

void test_ur3_fk_screw()
{
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}
    };

    std::cout << "UR3e FK using POE:\n\n";
    for (int i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d T = ur3e_fk_screw(joint_configurations[i]);
        print_pose("T", T);
    }
}

void test_ur3_fk_transform()
{
    std::vector<std::vector<double>> joint_configurations = {
        {0.0, 0.0, 0.0, -90.0, 0.0, 0.0},
        {0.0, -180.0, 0.0, 0.0, 0.0, 0.0},
        {0.0, -90.0, 0.0, 0.0, 0.0, 0.0}
    };

    std::cout << "UR3e FK using Homogenous Transformation:\n\n";
    for (int i = 0; i < joint_configurations.size(); i++) {
        Eigen::Matrix4d T = ur3e_fk_transform(joint_configurations[i]);
        print_pose("T", T);
    }
}

/*---------------- Here comes some test functions ----------------*/
// My own test functions since none where given
void test_Euler_zyx()
{
    constexpr double yaw = M_PI / 4;
    constexpr double pitch = - M_PI / 2;
    constexpr double roll = M_PI / 8;

    // Why did I make it take degrees input, in hindsight makes it really annoying to work with
    const Eigen::Vector3d euler_zyx_deg(rad_to_deg(yaw), rad_to_deg(pitch), rad_to_deg(roll));

    const Eigen::Matrix3d R = rotation_matrix_from_euler_zyx(euler_zyx_deg);

    const Eigen::Vector3d euler_zyx(yaw, pitch, roll);

    std::cout << "Original ZYX: " << euler_zyx.transpose() << std::endl;

    Eigen::Vector3d euler_zyx_recovered = euler_zyx_from_rotation_matrix(R);

    std::cout << "Recovered ZYX: " << euler_zyx_recovered.transpose() << "\n" << std::endl;
}

void test_twist_func()
{
    const Eigen::Vector3d w(10,20,30);
    const Eigen::Vector3d v(15,7,0);

    std::cout << "Original v: " << v.transpose() << " | Original w: " << w.transpose() << std::endl;
    std::cout << "New Twist Vector: " << twist(w, v).transpose() << "\n" << std::endl;
}

void test_screw_func()
{
    Eigen::Vector3d s(0, 0, 1);

    Eigen::Vector3d q(1, 0, 0);

    double h = 0.5;

    Eigen::VectorXd S = screw_axis(q, s, h);

    std::cout << "Screw axis S = " << S.transpose() << "\n" << std::endl;
}

void test_adjoint_matrix()
{
    Eigen::Matrix3d R = rotate_z(45);

    Eigen::Vector3d p (1,2,3);

    Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
    T.topLeftCorner<3,3>() = R;
    T.topRightCorner<3,1>() = p;

    Eigen::MatrixXd Adj = adjoint_matrix(T);

    std::cout << "AdjointMatrix: \n " << Adj << "\n" << std::endl;
}

// Could make test functions for further functions in task 3, however have not done it. I trust my code... famous last words

int main()
{
    /*test_Euler_zyx();
    test_twist_func();
    test_screw_func();
    test_adjoint_matrix(); */
    task2_a();
    task2_b();
    test_3r_planar_tf();
    test_3r_planar_screw();
    test_ur3_fk_screw();
    test_ur3_fk_transform();
    return 0;
}