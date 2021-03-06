#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <vector>


double NormaliseAngle(double angle);

namespace Lidar_index {
    constexpr int x = 0;
    constexpr int y = 1;
}

namespace Radar_index {
    constexpr int longitudinal_distance = 0;
    constexpr int angle_of_view = 1;
    constexpr int longitudinal_velocity = 2;
}

namespace UKF_index {
    constexpr int x = 0;
    constexpr int y = 1;
    constexpr int velocity = 2;
    constexpr int theta = 3;
    constexpr int theta_acc = 4;
    constexpr int mu_acc = 5;
    constexpr int mu_theta_acc_acc = 6;
}

class UKF {
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;
  std::vector<Eigen::VectorXd> x_history;
  Eigen::VectorXd x_aug_;

  // state covariance matrix
  Eigen::MatrixXd P_;
  Eigen::MatrixXd P_aug_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;
  Eigen::MatrixXd Xsig_pred_normalised;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  int n_sigma_;

  int n_z_radar_;
  int n_z_lidar_;

  // Sigma point spreading parameter
  double lambda_;

private:
    long long timestamp_;  // Most recently processed timestamp
    Eigen::MatrixXd sigma_points_;

    Eigen::VectorXd InitialiseFromRadar(MeasurementPackage meas_package);
    Eigen::VectorXd InitialiseFromLidar(MeasurementPackage meas_package);

    void EnsureCovarianceIsPositiveDefinite();
};

#endif  // UKF_H