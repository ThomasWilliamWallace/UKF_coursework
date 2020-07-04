#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

// Get rid of any negative / tiny covariance values
void UKF::EnsureCovarianceIsPositiveDefinite() {
//    for (int i = 0; i < P_.rows(); i++) {
//        for (int j = 0; j < P_.cols(); j++) {
//            if (P_(i, j) < 1e-6) {
//                P_(i, j) = 1e-6;
//            }
//        }
//    }
}

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
    // if this is false, laser measurements will be ignored (except during init)
    use_laser_ = true;

    // if this is false, radar measurements will be ignored (except during init)
    use_radar_ = true;

    n_x_ = 5;
    n_aug_ = 7;
    n_sigma_ = 1 + 2 * n_aug_;
    n_z_radar_ = 3;
    n_z_lidar_ = 2;

    // initial state vector
    x_ = VectorXd::Zero(n_x_);
    x_aug_ = VectorXd::Zero(n_aug_);

    // initial covariance matrix
    P_ = MatrixXd::Identity(n_x_, n_x_);
    P_(UKF_index::theta, UKF_index::theta) = 0.2;  // reduce sigma point spread so it doesn't wrap around a 2*PI interval. A smaller spread will have a better linearisation.
    P_aug_ = MatrixXd::Zero(n_aug_, n_aug_);

    Xsig_pred_ = MatrixXd::Zero(n_x_, n_sigma_);
    Xsig_pred_normalised = MatrixXd::Zero(n_x_, n_sigma_);

    // Process noise standard deviation longitudinal acceleration in m/s^2
    std_a_ = 3;

    // Process noise standard deviation yaw acceleration in rad/s^2
    std_yawdd_ = (M_PI * 3 / 6) / 2;

    /**
     * DO NOT MODIFY measurement noise values below.
     * These are provided by the sensor manufacturer.
     */

    // Laser measurement noise standard deviation position1 in m
    std_laspx_ = 0.15;

    // Laser measurement noise standard deviation position2 in m
    std_laspy_ = 0.15;

    // Radar measurement noise standard deviation radius in m
    std_radr_ = 0.3;

    // Radar measurement noise standard deviation angle in rad
    std_radphi_ = 0.03;

    // Radar measurement noise standard deviation radius change in m/s
    std_radrd_ = 0.3;

    /**
     * End DO NOT MODIFY section for measurement noise values
     */

    /**
     * TODO: Complete the initialization. See ukf.h for other member properties.
     * Hint: one or more values initialized above might be wildly off...
     */
    lambda_ = 3 - n_aug_;
    weights_ = Eigen::VectorXd(n_sigma_);
    weights_(0) = lambda_ / (lambda_ + n_aug_);
    for (int sigma = 1; sigma < weights_.size(); sigma++) {
        weights_(sigma) = 1 / (2 * (lambda_ + n_aug_));
    }
    sigma_points_ = Eigen::MatrixXd::Zero(n_aug_, n_sigma_);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
    /**
     * TODO: Complete this function! Make sure you switch between lidar and radar
     * measurements.
     */
    //  std::cout << "UKF::ProcessMeasurement, x_=" << x_ << "\n";
//    if (!is_initialized_ || abs(meas_package.timestamp_ - timestamp_) > 0) {
        Prediction((meas_package.timestamp_ - timestamp_) / 1000000.0);
//    }
    if (meas_package.sensor_type_ == meas_package.LASER) {
        UpdateLidar(meas_package);
    } else {
        UpdateRadar(meas_package);
    }
    timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
}

double NormaliseAngle(double angle) {
//    if (angle > 2 * M_PI) {
//        angle = std::fmod(angle, 2*M_PI);
//    }
    while (angle > M_PI) {
        angle -= 2 * M_PI;
    }
//    if (angle < 2 * M_PI) {
//        angle = -std::fmod(-angle, 2*M_PI);
//    }
    while (angle < -M_PI) {
        angle += 2 * M_PI;
    }
    return angle;
}

void UKF::Prediction(double delta_t) {
    /**
     * TODO: Complete this function! Estimate the object's location.
     * Modify the state vector, x_. Predict sigma points, the state,
     * and the state covariance matrix.
     */
//    std::cout << "\nUKF::Prediction\n";
    x_history.push_back(x_);
    if (!is_initialized_) {
        return;
    }
    std::cout << "Prediction Update\n";
//    std::cout << "x_=" << x_ << "\n";
//    std::cout << "weights_=\n" << weights_ << "\n";

    Eigen::VectorXd old_x = x_;
    x_aug_ << x_, 0, 0;
    P_aug_.block(0, 0, n_x_, n_x_) = P_;
    P_aug_(UKF_index::mu_acc, UKF_index::mu_acc) = std_a_ * std_a_;
    P_aug_(UKF_index::mu_theta_acc_acc, UKF_index::mu_theta_acc_acc) = std_yawdd_ * std_yawdd_;
    MatrixXd sqrt_matrix = P_aug_.llt().matrixL();
    MatrixXd sigma_diff = std::sqrt(lambda_ + n_aug_) * sqrt_matrix;
//    std::cout << "x_aug_=\n" << x_aug_ << "\n";
//    std::cout << "P_aug_=\n" << P_aug_ << "\n";
//    std::cout << "sqrt_matrix=\n" << sqrt_matrix << "\n";
//    std::cout << "sigma_diff=\n" << sigma_diff << "\n";

//    std::stringstream ss;
//    ss << x_aug_;
//    std::string str_x_aug_ = ss.str();
//    ss.str("");
//    ss << sigma_diff;
//    std::string str_sigma_diff = ss.str();
//    ss.str("");
//    std::string str_sigma_points;

    // Generate sigma points
    sigma_points_.fill(0.0);
    sigma_points_.col(0) = x_aug_;
//    ss << sigma_points_;
//    str_sigma_points = ss.str();
//    ss.str("");
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + sigma) = x_aug_ + sigma_diff.col(sigma);
//        ss << sigma_points_;
//        str_sigma_points = ss.str();
//        ss.str("");
        sigma_points_(UKF_index::theta, 1 + sigma) = NormaliseAngle(sigma_points_(UKF_index::theta, 1 + sigma));
//        ss << sigma_points_;
//        str_sigma_points = ss.str();
//        ss.str("");
    }
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + n_aug_ + sigma) = x_aug_ - sigma_diff.col(sigma);
//        ss << sigma_points_;
//        str_sigma_points = ss.str();
//        ss.str("");
        sigma_points_(UKF_index::theta, 1 + n_aug_ + sigma) = NormaliseAngle(sigma_points_(UKF_index::theta, 1 + n_aug_ + sigma));
//        ss << sigma_points_;
//        str_sigma_points = ss.str();
//        ss.str("");
    }
//    std::cout << "sigma_points_=\n" << sigma_points_ << "\n";
//    std::string str_P;
//    std::string str_x;
//    std::string str_Xsig_pred;
    Xsig_pred_.fill(0);
    double delta_t2 = delta_t * delta_t;
    // Apply process model to sigma points
    for (int sigma = 0; sigma < sigma_points_.cols(); sigma++) {
        bool angular_acceleration_is_zero;
        float x = sigma_points_(UKF_index::x, sigma);
        float y = sigma_points_(UKF_index::y, sigma);
        float velocity = sigma_points_(UKF_index::velocity, sigma);
        float theta = sigma_points_(UKF_index::theta, sigma);
        float theta_acc = sigma_points_(UKF_index::theta_acc, sigma);
        float mu_acc = sigma_points_(UKF_index::mu_acc, sigma);
        float mu_theta_acc_acc = sigma_points_(UKF_index::mu_theta_acc_acc, sigma);

        if (abs(theta_acc) < 1e-6) {
            angular_acceleration_is_zero = true;
        } else {
            angular_acceleration_is_zero = false;
        }
        if (angular_acceleration_is_zero) {
            Xsig_pred_(UKF_index::x, sigma) = x + delta_t * velocity * std::cos(theta) + delta_t2 * std::cos(theta) * mu_acc / 2.0;
        } else {
            Xsig_pred_(UKF_index::x, sigma) = x + (velocity / theta_acc) * (std::sin(theta + theta_acc * delta_t) - std::sin(theta)) + delta_t2 * std::cos(theta) * mu_acc / 2.0;
        }
//        ss << Xsig_pred_;
//        str_Xsig_pred = ss.str();
//        ss.str("");
        if (angular_acceleration_is_zero) {
            Xsig_pred_(UKF_index::y, sigma) = y + delta_t * velocity * std::sin(theta) + delta_t2 * std::sin(theta) * mu_acc / 2.0;
        } else {
            Xsig_pred_(UKF_index::y, sigma) = y + (velocity / theta_acc) * (-std::cos(theta + theta_acc * delta_t) + std::cos(theta)) + delta_t2 * std::sin(theta) * mu_acc / 2.0;
        }
//        ss << Xsig_pred_;
//        str_Xsig_pred = ss.str();
//        ss.str("");
        Xsig_pred_(UKF_index::velocity, sigma) = velocity + delta_t * mu_acc;
//        ss << Xsig_pred_;
//        str_Xsig_pred = ss.str();
//        ss.str("");
        Xsig_pred_(UKF_index::theta, sigma) = NormaliseAngle(theta + delta_t * theta_acc + delta_t2 * mu_theta_acc_acc / 2.0);
//        Xsig_pred_(UKF_index::theta, sigma) = theta + delta_t * theta_acc + delta_t2 * mu_theta_acc_acc / 2.0;
//        ss << Xsig_pred_;
//        str_Xsig_pred = ss.str();
//        ss.str("");
        Xsig_pred_(UKF_index::theta_acc, sigma) = theta_acc + delta_t * mu_theta_acc_acc;
//        ss << Xsig_pred_;
//        str_Xsig_pred = ss.str();
//        ss.str("");
    }
//    ss << P_;
//    str_P = ss.str();
//    ss.str("");
//    ss << x_;
//    str_x = ss.str();
//    ss.str("");
//    ss << Xsig_pred_;
//    str_Xsig_pred = ss.str();
//    ss.str("");
//    std::cout << "Xsig_pred_=\n" << Xsig_pred_ << "\n";

    // Calculate sigma point mean
    x_.fill(0);
//    ss << x_;
//    str_x = ss.str();
//    ss.str("");
    double mean_theta = Xsig_pred_(UKF_index::theta, 0   );
    // Centering the sigma centre point theta on 0, to allow averaging.
    Xsig_pred_normalised = Xsig_pred_;
    Xsig_pred_normalised.row(UKF_index::theta) -= mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
//    ss << Xsig_pred_normalised;
//    std::string str_Xsig_pred_normalised = ss.str();
//    ss.str("");
    for (int sigma = 0; sigma < n_sigma_; sigma++) {
        Xsig_pred_normalised(UKF_index::theta, sigma) = NormaliseAngle(Xsig_pred_normalised(UKF_index::theta, sigma));
    }
//    ss << Xsig_pred_normalised;
//    str_Xsig_pred_normalised = ss.str();
//    ss.str("");
    for (int state_index = 0; state_index < n_x_; state_index++) {
        for (int sig_index = 0; sig_index < n_sigma_; sig_index++) {
            x_(state_index) += Xsig_pred_normalised(state_index, sig_index) * weights_(sig_index);
//            ss << x_;
//            str_x = ss.str();
//            ss.str("");
        }
    }
    // Shifting the sigma point back to it's original value.
    Xsig_pred_normalised.row(UKF_index::theta) += mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
//    ss << Xsig_pred_;
//    str_Xsig_pred_normalised = ss.str();
//    ss.str("");
//    for (int sigma = 0; sigma < n_sigma_; sigma++) {
//        Xsig_pred_(UKF_index::theta, sigma) = NormaliseAngle(Xsig_pred_(UKF_index::theta, sigma));
//    }
//    ss << Xsig_pred_;
//    str_Xsig_pred = ss.str();
//    ss.str("");
    x_.row(UKF_index::theta) += mean_theta * Eigen::VectorXd::Ones(1).transpose();
//    ss << x_;
//    str_x = ss.str();
//    ss.str("");

    x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));
//    ss << x_;
//    str_x = ss.str();
//    ss.str("");
//    std::cout << "x_=\n" << x_ << "\n";

    // Calculate sigma point covariance
    P_.fill(0.0);
//    std::string str_Xsig_pred_col;
//    std::string str_weight;
//    std::string str_diff;
//    std::string str_diff_squared;
//    std::string str_contribution;
    for (int sig_index = 0; sig_index < n_sigma_; sig_index++) {
//        ss << P_;
//        str_P = ss.str();
//        ss.str("");
//        ss << x_;
//        str_x = ss.str();
//        ss.str("");
//        ss << Xsig_pred_normalised.col(sig_index);
//        str_Xsig_pred_col = ss.str();
//        ss.str("");
//        std::stringstream ss_weight;
//        ss_weight << weights_(sig_index);
//        str_weight = ss_weight.str();
        VectorXd difference = Xsig_pred_.col(sig_index) - x_;
//        std::stringstream ss_diff;
//        ss_diff << difference;
//        str_diff = ss_diff.str();
        MatrixXd diff_squared = difference * difference.transpose();
//        std::stringstream ss_diff_squared;
//        ss_diff_squared << diff_squared;
//        str_diff_squared = ss_diff_squared.str();
        MatrixXd contribution = diff_squared * weights_(sig_index);
//        std::stringstream ss_contribution;
//        ss_contribution << contribution;
//        str_contribution = ss_contribution.str();
        P_ += difference * difference.transpose() * weights_(sig_index);
//        ss.str("");
//        ss << P_;
//        str_P = ss.str();
//        ss.str("");
    }

    EnsureCovarianceIsPositiveDefinite();

    std::cout << "x_=\n" << x_ << "\n";
    std::cout << "P_=\n" << P_ << "\n";
    double delta_theta = x_(UKF_index::theta) - old_x(UKF_index::theta);
    delta_theta += (delta_theta>M_PI) ? -M_PI*2 : (delta_theta<-M_PI) ? M_PI*2 : 0;
//    if (abs(delta_theta) > 1.6) {
//        ss << old_x;
//        std::string str_old_x = ss.str();
//        ss.str("");
//    }
}

Eigen::VectorXd UKF::LidarMeasurementFunction(MeasurementPackage meas_package) {
    double meas_x = meas_package.raw_measurements_(Lidar_index::x);
    double meas_y = meas_package.raw_measurements_(Lidar_index::y);

    double meas_forward_velocity = 0; // NOT MEASURED
    double meas_theta = 0;  // NOT MEASURED
    double meas_ang_acc = 0; // NOT MEASURED

    Eigen::VectorXd measured_state = Eigen::VectorXd(n_x_);
    measured_state << meas_x, meas_y, meas_forward_velocity, meas_theta, meas_ang_acc;
    return measured_state;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
    /**
     * TODO: Complete this function! Use lidar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the lidar NIS, if desired.
     */
    Eigen::VectorXd measured_state = LidarMeasurementFunction(meas_package);

    if (is_initialized_) {
        std::cout << "Lidar Measurement Update\n";
//        std::stringstream ss;
        Eigen::VectorXd old_x = x_;
        // Apply measurement model to sigma points
        Eigen::MatrixXd z_sig = Eigen::MatrixXd(n_z_lidar_, n_sigma_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            double x = Xsig_pred_(UKF_index::x, sigma);
            double y = Xsig_pred_(UKF_index::y, sigma);
            z_sig(Lidar_index::x, sigma) = x;
            z_sig(Lidar_index::y, sigma) = y;
        }
//        ss << z_sig;
//        std::string str_z_sig = ss.str();
//        ss.str("");

        // Calculate mean and covariance from sigma points
        Eigen::VectorXd z_mean = Eigen::VectorXd::Zero(n_z_lidar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            for (int row = 0; row < z_mean.size(); row++) {
                z_mean(row) += z_sig(row, sigma) * weights_(sigma);
            }
        }
//        ss << z_mean;
//        std::string str_z_mean = ss.str();
//        ss.str("");

        MatrixXd R = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
        R(Lidar_index::x, Lidar_index::x) = std_laspx_ * std_laspx_;
        R(Lidar_index::y, Lidar_index::y) = std_laspy_ * std_laspy_;
//        ss << R;
//        std::string str_R = ss.str();
//        ss.str("");

        MatrixXd S = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
        S += R;
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
//            std::cout << "sigma=" << sigma << "\n";
            Eigen::VectorXd difference = z_sig.col(sigma) - z_mean;
//            ss << difference;
//            std::string str_z_sig_z_mean_diff = ss.str();
//            std::cout << "str_z_sig_z_mean_diff=\n" << str_z_sig_z_mean_diff << "\n";
//            ss.str("");
            S += difference * difference.transpose() * weights_(sigma);
        }
//        ss << S;
//        std::string str_S = ss.str();
//        ss.str("");

        // UKF update (can be shared with radar)
        Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_x_, n_z_lidar_);
//        ss << Tc;
//        std::string str_Tc = ss.str();
//        ss.str("");
//        ss << Xsig_pred_;
//        std::string str_Xsig_pred_ = ss.str();
//        ss.str("");
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
//            std::cout << "sigma=" << sigma << "\n";
            VectorXd x_diff = Xsig_pred_normalised.col(sigma) - x_;
//            ss << x_diff;
//            std::string str_tc_x_diff = ss.str();
//            std::cout << "str_tc_x_diff=\n" << str_tc_x_diff << "\n";
//            ss.str("");
            x_diff(UKF_index::theta) = NormaliseAngle(x_diff(UKF_index::theta));
//            ss << x_diff;
//            std::string str_tc_x_diff_normalised = ss.str();
//            std::cout << "tc_x_diff_normalised=\n" << str_tc_x_diff_normalised << "\n";
//            ss.str("");
            VectorXd z_diff = z_sig.col(sigma) - z_mean;
//            ss << z_diff;
//            std::string str_tc_z_diff = ss.str();
//            std::cout << "str_tc_z_diff=\n" << str_tc_z_diff << "\n";
//            ss.str("");
            Eigen::MatrixXd mult = x_diff * z_diff.transpose();
//            ss << mult;
//            std::string str_tc_mult = ss.str();
//            std::cout << "str_tc_mult=\n";
//            std::cout << str_tc_mult << "\n";
//            ss.str("");
            Tc = Tc + weights_(sigma) * mult;
//            ss << Tc;
//            str_Tc = ss.str();
//            ss.str("");
//            std::cout << "\n";
        }

        // calculate Kalman gain K;
//        ss << S.inverse();
//        std::string str_S_inverse = ss.str();
//        ss.str("");
        MatrixXd K;
        K = Tc * S.inverse();
//        ss << K;
//        std::string str_K = ss.str();
//        ss.str("");

        // update state mean and covariance matrix
//        ss << meas_package.raw_measurements_;
//        std::string str_meas = ss.str();
//        ss.str("");
        VectorXd z_diff = meas_package.raw_measurements_ - z_mean;
//        ss << z_diff;
//        std::string str_z_diff = ss.str();
//        ss.str("");

//        ss << K * z_diff;
//        std::string str_k_z_diff = ss.str();
//        ss.str("");
        x_ = x_ + K * z_diff;
//        ss << x_;
//        std::string str_x = ss.str();
//        ss.str("");
        x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));
//        ss << x_;
//        std::string str_x_normalised = ss.str();
//        ss.str("");
        double delta_theta = x_(UKF_index::theta) - old_x(UKF_index::theta);
        delta_theta += (delta_theta>M_PI) ? -M_PI*2 : (delta_theta<-M_PI) ? M_PI*2 : 0;
//        if (abs(delta_theta) > 1.6) {
//            ss << old_x;
//            std::string str_old_x = ss.str();
//            ss.str("");
//        }
        P_ = P_ - K * S * K.transpose();

        EnsureCovarianceIsPositiveDefinite();

        std::cout << "x_=\n" << x_ << "\n";
        std::cout << "P_=\n" << P_ << "\n";

    } else {
        x_ << measured_state;
        x_(UKF_index::velocity) = 5;
//        x_(UKF_index::theta) = M_PI/2;
//      x_(UKF_index::theta_acc) = M_PI/5;
//        P_ *= 1e-6;
    }
}

Eigen::VectorXd UKF::RadarMeasurementFunction(MeasurementPackage meas_package) {
    double meas_longitudinal_distance = meas_package.raw_measurements_(Radar_index::longitudinal_distance);
    double meas_positional_theta = meas_package.raw_measurements_(Radar_index::angle_of_view);
    double meas_longitudinal_velocity = meas_package.raw_measurements_(Radar_index::longitudinal_velocity);

    double meas_x = std::cos(meas_positional_theta) * meas_longitudinal_distance;
    double meas_y = std::sin(meas_positional_theta) * meas_longitudinal_distance;
    double meas_forward_velocity = 0;  // NOT MEASURED
    double meas_theta = 0;  // NOT MEASURED
    double meas_ang_acc = 0;  // NOT MEASURED

    Eigen::VectorXd measured_state = Eigen::VectorXd(n_x_);
    measured_state << meas_x, meas_y, meas_forward_velocity, meas_theta, meas_ang_acc;
    return measured_state;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
    /**
     * TODO: Complete this function! Use radar data to update the belief
     * about the object's position. Modify the state vector, x_, and
     * covariance, P_.
     * You can also calculate the radar NIS, if desired.
     */
    Eigen::VectorXd measured_state = RadarMeasurementFunction(meas_package);

    if (is_initialized_) {
        std::cout << "Radar Measurement Update\n";
//        std::stringstream ss;
//        ss << x_;
//        std::string str_x_ = ss.str();
//        ss.str("");
//        ss << Xsig_pred_;
//        std::string str_Xsig_pred_ = ss.str();
//        ss.str("");
//        ss << Xsig_pred_normalised;
//        std::string str_Xsig_pred_normalised = ss.str();
//        ss.str("");
        Eigen::VectorXd old_x = x_;
        // Apply measurement model to sigma points
        Eigen::MatrixXd z_sig = Eigen::MatrixXd(n_z_radar_, n_sigma_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            double x = Xsig_pred_normalised(UKF_index::x, sigma);
            double y = Xsig_pred_normalised(UKF_index::y, sigma);
            double velocity = Xsig_pred_normalised(UKF_index::velocity, sigma);
            double theta = Xsig_pred_normalised(UKF_index::theta, sigma);
            double longitudinal_distance = std::sqrt(x * x + y * y);
            double angle_of_view = std::atan2(y, x);
            double longitudinal_velocity = (x * std::cos(theta) * velocity + y * std::sin(theta) * velocity) / longitudinal_distance;
            z_sig(Radar_index::longitudinal_distance, sigma) = longitudinal_distance;
            z_sig(Radar_index::angle_of_view, sigma) = angle_of_view;
            z_sig(Radar_index::longitudinal_velocity, sigma) = longitudinal_velocity;
        }
//        ss << z_sig;
//        std::string str_z_sig = ss.str();
//        ss.str("");

        // Calculate mean and covariance from sigma points
        // center angle on the sigma mean
        double mean_theta = z_sig(Radar_index::angle_of_view, 0  );
        z_sig.row(Radar_index::angle_of_view) -= mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            z_sig(Radar_index::angle_of_view, sigma) = NormaliseAngle(z_sig(Radar_index::angle_of_view, sigma));
        }
        Eigen::VectorXd z_mean = Eigen::VectorXd::Zero(n_z_radar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            for (int row = 0; row < z_mean.size(); row++) {
                z_mean(row) += z_sig(row, sigma) * weights_(sigma);
            }
        }
        z_mean(Radar_index::angle_of_view) = NormaliseAngle(z_mean(1));
        // uncenter angle from the sigma mean
        z_sig.row(Radar_index::angle_of_view) += mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
        z_mean(Radar_index::angle_of_view) += mean_theta;
        z_mean(Radar_index::angle_of_view) = NormaliseAngle(z_mean(Radar_index::angle_of_view));
//        ss << z_mean;
//        std::string str_z_mean = ss.str();
//        ss.str("");

        MatrixXd R = MatrixXd::Zero(n_z_radar_, n_z_radar_);
        R(Radar_index::longitudinal_distance, Radar_index::longitudinal_distance) = std_radr_ * std_radr_;
        R(Radar_index::angle_of_view, Radar_index::angle_of_view) = std_radphi_ * std_radphi_;
        R(Radar_index::longitudinal_velocity, Radar_index::longitudinal_velocity) = std_radrd_ * std_radrd_;
//        ss << R;
//        std::string str_R = ss.str();
//        ss.str("");

        MatrixXd S = MatrixXd::Zero(n_z_radar_, n_z_radar_);
        S += R;
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            Eigen::VectorXd difference = z_sig.col(sigma) - z_mean;
            S += difference * difference.transpose() * weights_(sigma);
        }
//        ss << S;
//        std::string str_S = ss.str();
//        ss.str("");

        // UKF update (can be shared with lidar)
        Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_x_, n_z_radar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            VectorXd x_diff = Xsig_pred_normalised.col(sigma) - x_;
            x_diff(UKF_index::theta) = NormaliseAngle(x_diff(UKF_index::theta));
            VectorXd z_diff = z_sig.col(sigma) - z_mean;
            z_diff(Radar_index::angle_of_view) = NormaliseAngle(z_diff(Radar_index::angle_of_view));
            Tc = Tc + weights_(sigma) * x_diff * z_diff.transpose();
        }
//        ss << Tc;
//        std::string str_Tc = ss.str();
//        ss.str("");

        // calculate Kalman gain K;
//        ss << S.inverse();
//        std::string str_S_inverse = ss.str();
//        ss.str("");
        MatrixXd K;
        K = Tc * S.inverse();
//        ss << K;
//        std::string str_K = ss.str();
//        ss.str("");

        // update state mean and covariance matrix
//        ss << meas_package.raw_measurements_;
//        std::string str_meas = ss.str();
//        ss.str("");
        VectorXd z_diff = meas_package.raw_measurements_ - z_mean;
//        ss << z_diff;
//        std::string str_z_diff = ss.str();
//        ss.str("");
        z_diff(Radar_index::angle_of_view) = NormaliseAngle(z_diff(Radar_index::angle_of_view));
//        ss << z_diff;
//        std::string str_z_diff_normalised = ss.str();
//        ss.str("");

//        ss << K * z_diff;
//        std::string str_k_z_diff = ss.str();
//        ss.str("");
        x_ = x_ + K * z_diff;
//        ss << x_;
//        std::string str_x = ss.str();
//        ss.str("");
        x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));
//        ss << x_;
//        std::string str_x_normalised = ss.str();
//        ss.str("");
        double delta_theta = x_(UKF_index::theta) - old_x(UKF_index::theta);
        delta_theta += (delta_theta>M_PI) ? -M_PI*2 : (delta_theta<-M_PI) ? M_PI*2 : 0;
        if (abs(delta_theta) > 1.6) {
//            ss << old_x;
//            std::string str_old_x = ss.str();
//            ss.str("");
        }
        P_ = P_ - K * S * K.transpose();

        EnsureCovarianceIsPositiveDefinite();

        std::cout << "x_=\n" << x_ << "\n";
        std::cout << "P_=\n" << P_ << "\n";

    } else {
        x_ << measured_state;
        x_(UKF_index::velocity) = 5;
//        x_(UKF_index::theta) = M_PI/2;
//      x_(UKF_index::theta_acc) = M_PI/5;
//      P_ << 0;
//        P_ *= 1e-6;
    }
}