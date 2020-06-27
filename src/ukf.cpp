#include "ukf.h"
#include "Eigen/Dense"
#include <cmath>
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // initial state vector
  x_ = VectorXd(n_x_);
  x_aug_ = VectorXd(n_aug_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_aug_ = MatrixXd(n_aug_, n_aug_);
  P_aug_.fill(0.0);

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
   for (int sigma = 0; sigma < weights_.size(); sigma++) {
       weights_(sigma) = 1 / (2 * (lambda_ + n_aug_));
   }
   sigma_points_ = Eigen::MatrixXd(n_aug_, n_sigma_);
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  //  std::cout << "UKF::ProcessMeasurement, x_=" << x_ << "\n";
  Prediction(meas_package.timestamp_ - timestamp_);
  if (meas_package.sensor_type_ == meas_package.LASER) {
      UpdateLidar(meas_package);
  } else {
      UpdateRadar(meas_package);
  }
  timestamp_ = meas_package.timestamp_;
  is_initialized_ = true;
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
   std::cout << "UKF::Prediction, x_=" << x_ << "\n";

   x_aug_ << x_, 0, 0;
   P_aug_.block(0, 0, n_x_, n_x_) = P_;
   P_aug_(n_x_, n_x_) = std_a_ * std_a_;
   P_aug_(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;
   MatrixXd sqrt_matrix = P_aug_.llt().matrixL();
   MatrixXd sigma_diff = std::sqrt(lambda_ + n_aug_) * sqrt_matrix;

    // Generate sigma points
    sigma_points_.col(0) = x_aug_;
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + sigma) = x_aug_ + sigma_diff.col(sigma);
    }
    for (int sigma = 0; sigma < n_aug_; sigma++) {
        sigma_points_.col(1 + n_aug_ + sigma) = x_aug_ - sigma_diff.col(sigma);
    }

    // Apply process model to sigma points


    // Calculate mean and covariance of sigma points
}

Eigen::VectorXd UKF::LidarMeasurementFunction(MeasurementPackage meas_package) {
    double meas_x = meas_package.raw_measurements_(0);
    double meas_y = meas_package.raw_measurements_(1);

    double meas_forward_velocity = 0; // NOT MEASURED
    double meas_theta = 0;  // NOT MEASURED
    double meas_ang_acc = 0; // NOT MEASURED

    Eigen::VectorXd measured_state = Eigen::VectorXd(5);
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
        // Apply measurement model to sigma points

        // Calculate mean and covariance from sigma points

    } else {
        x_ << measured_state;
//        P_ << ;
    }
}

Eigen::VectorXd UKF::RadarMeasurementFunction(MeasurementPackage meas_package) {
    double meas_longitudinal_distance = meas_package.raw_measurements_(0);
    double meas_positional_theta = meas_package.raw_measurements_(1);
    double meas_longitudinal_velocity = meas_package.raw_measurements_(2);

    double meas_x = std::cos(meas_positional_theta) * meas_longitudinal_distance;
    double meas_y = std::sin(meas_positional_theta) * meas_longitudinal_distance;
    double meas_forward_velocity = 0;  // NOT MEASURED
    double meas_theta = 0;  // NOT MEASURED
    double meas_ang_acc = 0;  // NOT MEASURED

    Eigen::VectorXd measured_state = Eigen::VectorXd(5);
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
      // Apply measurement model to sigma points

      // Calculate mean and covariance from sigma points

  } else {
      x_ << measured_state;
//      P_ << 0;
  }
}