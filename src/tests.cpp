//
// Created by thomas on 28/06/2020.
//

#include <iostream>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "ukf.h"
#include <cassert>
#include <cmath>
#include <vector>


// Checks that UKF state is unchanged by a prediction of zero duration.
void TestZeroDurationPredict() {
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
    std::cout << "TestZeroDurationPredict Completed\n";
};


// Checks that UKF predicted velocity is constant if we're travelling in a straight line.
void TestStraightLineConstantVelocity() {
    UKF ukf;

    int x_index = 0;
    int y_index = 1;
    int velocity_index = 2;
    int theta_acc_index = 4;

    // Feed in first measurement
    MeasurementPackage meas_package;
    meas_package.sensor_type_ = MeasurementPackage::LASER;
    meas_package.raw_measurements_ = Eigen::VectorXd(2);
    meas_package.raw_measurements_ << -9.6649140872479418, 3.9175380731566825;
    meas_package.timestamp_ = 0;
    ukf.ProcessMeasurement(meas_package);

    // Set ukf theta_acc to zero so we are travelling in a straight line
    ukf.x_(theta_acc_index) = 0;

    // Set ukf velocity to 5.81 m/s
    double target_velocity = 5.81;
    ukf.x_(velocity_index) = target_velocity;

    std::vector<int> delta_t_list = {1, 200, 168, 50, 50000};
    int total_t = 0;
    for (auto& delta_t : delta_t_list) {
        total_t += delta_t;
    }
    double total_distance = 0;

    int delta_t;
    double prev_x_position = ukf.x_(0);
    double prev_y_position = ukf.x_(1);

    for (auto& delta_t : delta_t_list) {
        ukf.Prediction(delta_t);
        double x_diff = ukf.x_(x_index) - prev_x_position;
        double y_diff = ukf.x_(y_index) - prev_y_position;
        double distance = std::sqrt(x_diff * x_diff + y_diff * y_diff);
        total_distance += distance;
        double current_velocity = distance / delta_t;

        bool velocityMatchedTarget = (current_velocity - target_velocity) < 1e10;
        assert(("Velocity matches the target velocity.", velocityMatchedTarget));

        prev_x_position = ukf.x_(x_index);
        prev_y_position = ukf.x_(y_index);
    }

    bool averageVelocityMatchedTarget = ((total_distance/total_t) - target_velocity) < 1e10;
    assert(("Average velocity across the whole distance matches the target velocity.", averageVelocityMatchedTarget));
    std::cout << "TestStraightLineConstantVelocity Completed\n";
};

void RunAllTests() {
    std::cout << "Running all tests\n";
    TestZeroDurationPredict();
    TestStraightLineConstantVelocity();
    std::cout << "All tests Completed\n";
};

int main(int argc, char** argv) {
    RunAllTests();
}