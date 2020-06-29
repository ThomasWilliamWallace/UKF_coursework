//
// Created by thomas on 28/06/2020.
//

#include <iostream>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "ukf.h"
#include <cassert>


// Checks that UKF state is unchanged by a prediction of zero duration.
bool ZeroDurationPredict() {
    UKF ukf;

    // Feed in first measurement
    MeasurementPackage meas_package;
    meas_package.sensor_type_ = MeasurementPackage::LASER;
    meas_package.raw_measurements_ = Eigen::VectorXd(2);
    meas_package.raw_measurements_ << -9.6649140872479418, 3.9175380731566825;
    meas_package.timestamp_ = 0;
    ukf.ProcessMeasurement(meas_package);

    // Store UKF state
    Eigen::VectorXd x_prev = ukf.x_;
    Eigen::MatrixXd P_prev = ukf.P_;

    // Predict state after 0 time
    ukf.Prediction(0);

    Eigen::VectorXd x_diff = x_prev - ukf.x_;
    Eigen::MatrixXd P_diff = P_prev - ukf.P_;

    bool x_mean_changed = x_diff.isMuchSmallerThan(1e7);
    assert(("Mean state unchanged between initialization and zero-duration predict step.", x_mean_changed));

    bool P_mean_changed = P_diff.isMuchSmallerThan(1e7);
    assert(("Mean covariance matrix unchanged between initialization and zero-duration predict step.", P_mean_changed));

    // Store UKF state
    x_prev = ukf.x_;
    P_prev = ukf.P_;

    // Predict state after 0 time
    ukf.Prediction(0);

    x_diff = x_prev - ukf.x_;
    P_diff = P_prev - ukf.P_;

    x_mean_changed = x_diff.isMuchSmallerThan(1e10);
    assert(("Mean state unchanged between first and second zero-duration predict step.", x_mean_changed));

    P_mean_changed = P_diff.isMuchSmallerThan(1e10);
    assert(("Mean covariance matrix unchanged between first and second zero-duration predict step.", P_mean_changed));

    return false;
};

void RunAllTests() {
    std::cout << "Running all tests\n";
    ZeroDurationPredict();
};

int main(int argc, char** argv) {
    RunAllTests();
}