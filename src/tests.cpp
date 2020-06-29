//
// Created by thomas on 28/06/2020.
//

#include <iostream>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "ukf.h"


// Checks that UKF state is unchanged by a prediction of zero duration.
bool ZeroTimePredict() {
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

    // Verify that the UKF state is unchanged
    Eigen::VectorXd x_diff = x_prev - ukf.x_;
    Eigen::MatrixXd P_diff = P_prev - ukf.P_;
    std::cout << "x_prev=\n" << x_prev << "\n";
    std::cout << "ukf.x_=\n" << ukf.x_ << "\n";
    std::cout << "x_diff=\n" << x_diff << "\n";
    std::cout << "P_prev=\n" << P_prev << "\n";
    std::cout << "ukf.P_=\n" << ukf.P_ << "\n";
    std::cout << "P_diff=\n" << P_diff << "\n";

    return false;
};

void RunAllTests() {
    std::cout << "Running all tests\n";
    ZeroTimePredict();
};

int main(int argc, char** argv) {
    RunAllTests();
}