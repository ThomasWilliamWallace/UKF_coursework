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
//    P_ *= 100;
    P_(UKF_index::x, UKF_index::x) = 1;
    P_(UKF_index::y, UKF_index::y) = 1;
    P_(UKF_index::velocity, UKF_index::velocity) = 10;
    P_(UKF_index::theta, UKF_index::theta) = 0.2;  // reduce sigma point spread so it doesn't wrap around a 2*PI interval. A smaller spread will have a better linearisation.
//    P_(UKF_index::theta_acc, UKF_index::theta_acc) = 3;
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
        Prediction((meas_package.timestamp_ - timestamp_) / 1000000.0);
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

/**
 * Modify the state vector, x_. Predict sigma points, the state,
 * and the state covariance matrix.
 */
void UKF::Prediction(double delta_t) {
    x_history.push_back(x_);
    if (!is_initialized_) {
        return;
    }

    Eigen::VectorXd old_x = x_;
    x_aug_ << x_, 0, 0;
    P_aug_.block(0, 0, n_x_, n_x_) = P_;
    P_aug_(UKF_index::mu_acc, UKF_index::mu_acc) = std_a_ * std_a_;
    P_aug_(UKF_index::mu_theta_acc_acc, UKF_index::mu_theta_acc_acc) = std_yawdd_ * std_yawdd_;
    MatrixXd sqrt_matrix = P_aug_.llt().matrixL();
    MatrixXd sigma_diff = std::sqrt(lambda_ + n_aug_) * sqrt_matrix;

    // Generate sigma points
    sigma_points_.fill(0.0);
    sigma_points_.col(0) = x_aug_;
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + sigma) = x_aug_ + sigma_diff.col(sigma);
        sigma_points_(UKF_index::theta, 1 + sigma) = NormaliseAngle(sigma_points_(UKF_index::theta, 1 + sigma));
    }
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + n_aug_ + sigma) = x_aug_ - sigma_diff.col(sigma);
        sigma_points_(UKF_index::theta, 1 + n_aug_ + sigma) = NormaliseAngle(sigma_points_(UKF_index::theta, 1 + n_aug_ + sigma));
    }
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
        if (angular_acceleration_is_zero) {
            Xsig_pred_(UKF_index::y, sigma) = y + delta_t * velocity * std::sin(theta) + delta_t2 * std::sin(theta) * mu_acc / 2.0;
        } else {
            Xsig_pred_(UKF_index::y, sigma) = y + (velocity / theta_acc) * (-std::cos(theta + theta_acc * delta_t) + std::cos(theta)) + delta_t2 * std::sin(theta) * mu_acc / 2.0;
        }
        Xsig_pred_(UKF_index::velocity, sigma) = velocity + delta_t * mu_acc;
        Xsig_pred_(UKF_index::theta, sigma) = NormaliseAngle(theta + delta_t * theta_acc + delta_t2 * mu_theta_acc_acc / 2.0);
        Xsig_pred_(UKF_index::theta_acc, sigma) = theta_acc + delta_t * mu_theta_acc_acc;
    }

    // Calculate sigma point mean
    x_.fill(0);
    double mean_theta = Xsig_pred_(UKF_index::theta, 0   );
    // Centering the sigma centre point theta on 0, to allow averaging.
    Xsig_pred_normalised = Xsig_pred_;
    Xsig_pred_normalised.row(UKF_index::theta) -= mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
    for (int sigma = 0; sigma < n_sigma_; sigma++) {
        Xsig_pred_normalised(UKF_index::theta, sigma) = NormaliseAngle(Xsig_pred_normalised(UKF_index::theta, sigma));
    }
    for (int state_index = 0; state_index < n_x_; state_index++) {
        for (int sig_index = 0; sig_index < n_sigma_; sig_index++) {
            x_(state_index) += Xsig_pred_normalised(state_index, sig_index) * weights_(sig_index);
        }
    }
    // Shifting the sigma point back to it's original value.
    Xsig_pred_normalised.row(UKF_index::theta) += mean_theta * Eigen::VectorXd::Ones(n_sigma_).transpose();
    x_.row(UKF_index::theta) += mean_theta * Eigen::VectorXd::Ones(1).transpose();

    x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));

    // Calculate sigma point covariance
    P_.fill(0.0);
    for (int sig_index = 0; sig_index < n_sigma_; sig_index++) {
        VectorXd difference = Xsig_pred_.col(sig_index) - x_;
        MatrixXd diff_squared = difference * difference.transpose();
        MatrixXd contribution = diff_squared * weights_(sig_index);
        P_ += difference * difference.transpose() * weights_(sig_index);
    }

    EnsureCovarianceIsPositiveDefinite();

}

Eigen::VectorXd UKF::InitialiseFromLidar(MeasurementPackage meas_package) {
    double meas_x = meas_package.raw_measurements_(Lidar_index::x);
    double meas_y = meas_package.raw_measurements_(Lidar_index::y);

    double meas_forward_velocity = 0; // NOT MEASURED
    double meas_theta = 0;  // NOT MEASURED
    double meas_ang_acc = 0; // NOT MEASURED

    Eigen::VectorXd measured_state = Eigen::VectorXd(n_x_);
    measured_state << meas_x, meas_y, meas_forward_velocity, meas_theta, meas_ang_acc;
    return measured_state;
}

/**
 * Use lidar data to update the belief about the object's position.
 * Modify the state vector, x_, and covariance, P_.
 * You can also calculate the lidar NIS, if desired.
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
    if (is_initialized_) {
        Eigen::VectorXd old_x = x_;
        // Apply measurement model to sigma points
        Eigen::MatrixXd z_sig = Eigen::MatrixXd(n_z_lidar_, n_sigma_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            double x = Xsig_pred_(UKF_index::x, sigma);
            double y = Xsig_pred_(UKF_index::y, sigma);
            z_sig(Lidar_index::x, sigma) = x;
            z_sig(Lidar_index::y, sigma) = y;
        }

        // Calculate mean and covariance from sigma points
        Eigen::VectorXd z_mean = Eigen::VectorXd::Zero(n_z_lidar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            for (int row = 0; row < z_mean.size(); row++) {
                z_mean(row) += z_sig(row, sigma) * weights_(sigma);
            }
        }

        MatrixXd R = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
        R(Lidar_index::x, Lidar_index::x) = std_laspx_ * std_laspx_;
        R(Lidar_index::y, Lidar_index::y) = std_laspy_ * std_laspy_;

        MatrixXd S = MatrixXd::Zero(n_z_lidar_, n_z_lidar_);
        S += R;
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            Eigen::VectorXd difference = z_sig.col(sigma) - z_mean;
            S += difference * difference.transpose() * weights_(sigma);
        }

        // UKF update
        Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_x_, n_z_lidar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            VectorXd x_diff = Xsig_pred_normalised.col(sigma) - x_;
            x_diff(UKF_index::theta) = NormaliseAngle(x_diff(UKF_index::theta));
            VectorXd z_diff = z_sig.col(sigma) - z_mean;
            Eigen::MatrixXd mult = x_diff * z_diff.transpose();
            Tc = Tc + weights_(sigma) * mult;
        }

        // calculate Kalman gain K;
        MatrixXd K;
        K = Tc * S.inverse();

        // update state mean and covariance matrix
        VectorXd z_diff = meas_package.raw_measurements_ - z_mean;
        x_ = x_ + K * z_diff;
        x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));
        P_ = P_ - K * S * K.transpose();

        EnsureCovarianceIsPositiveDefinite();

    } else {
        x_ << InitialiseFromLidar(meas_package);
    }
}

Eigen::VectorXd UKF::InitialiseFromRadar(MeasurementPackage meas_package) {
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

/**
 * Use radar data to update the belief
 * about the object's position. Modify the state vector, x_, and
 * covariance, P_.
 * You can also calculate the radar NIS, if desired.
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

    if (is_initialized_) {
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

        MatrixXd R = MatrixXd::Zero(n_z_radar_, n_z_radar_);
        R(Radar_index::longitudinal_distance, Radar_index::longitudinal_distance) = std_radr_ * std_radr_;
        R(Radar_index::angle_of_view, Radar_index::angle_of_view) = std_radphi_ * std_radphi_;
        R(Radar_index::longitudinal_velocity, Radar_index::longitudinal_velocity) = std_radrd_ * std_radrd_;

        MatrixXd S = MatrixXd::Zero(n_z_radar_, n_z_radar_);
        S += R;
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            Eigen::VectorXd difference = z_sig.col(sigma) - z_mean;
            S += difference * difference.transpose() * weights_(sigma);
        }

        // UKF update (can be shared with lidar)
        Eigen::MatrixXd Tc = Eigen::MatrixXd::Zero(n_x_, n_z_radar_);
        for (int sigma = 0; sigma < n_sigma_; sigma++) {
            VectorXd x_diff = Xsig_pred_normalised.col(sigma) - x_;
            x_diff(UKF_index::theta) = NormaliseAngle(x_diff(UKF_index::theta));
            VectorXd z_diff = z_sig.col(sigma) - z_mean;
            z_diff(Radar_index::angle_of_view) = NormaliseAngle(z_diff(Radar_index::angle_of_view));
            Tc = Tc + weights_(sigma) * x_diff * z_diff.transpose();
        }

        // calculate Kalman gain K;
        MatrixXd K;
        K = Tc * S.inverse();

        // update state mean and covariance matrix
        VectorXd z_diff = meas_package.raw_measurements_ - z_mean;
        z_diff(Radar_index::angle_of_view) = NormaliseAngle(z_diff(Radar_index::angle_of_view));
        x_ = x_ + K * z_diff;
        x_(UKF_index::theta) = NormaliseAngle(x_(UKF_index::theta));
        P_ = P_ - K * S * K.transpose();

        EnsureCovarianceIsPositiveDefinite();

    } else {
        x_ << InitialiseFromRadar(meas_package);
    }
}