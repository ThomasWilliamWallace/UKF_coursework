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

    // Feed in first measurement
    MeasurementPackage meas_package;
    meas_package.sensor_type_ = MeasurementPackage::LASER;
    meas_package.raw_measurements_ = Eigen::VectorXd(2);
    meas_package.raw_measurements_ << -9.6649140872479418, 3.9175380731566825;
    meas_package.timestamp_ = 0;
    ukf.ProcessMeasurement(meas_package);

    // Set ukf theta_acc to zero so we are travelling in a straight line
    ukf.x_(UKF_index::theta_acc) = 0;

    // Set ukf velocity to 5.81 m/s
    double target_velocity = 5.81;
    ukf.x_(UKF_index::velocity) = target_velocity;

    std::vector<double> delta_t_list = {1, 168, 50000};
    for (int i = 0; i < delta_t_list.size(); i++) {
        delta_t_list[i] = delta_t_list[i] / 1000000.0;
    }
    int total_t = 0;
    for (auto& delta_t : delta_t_list) {
        total_t += delta_t;
    }
    double total_distance = 0;

    double prev_x_position = ukf.x_(UKF_index::x);
    double prev_y_position = ukf.x_(UKF_index::y);

    for (auto& delta_t : delta_t_list) {
        ukf.Prediction(std::round(delta_t * 1000000));
        double x_diff = ukf.x_(UKF_index::x) - prev_x_position;
        double y_diff = ukf.x_(UKF_index::y) - prev_y_position;
        double distance = std::sqrt(x_diff * x_diff + y_diff * y_diff);
        total_distance += distance;
        double current_velocity = distance / delta_t;

        bool velocityMatchedTarget = abs(current_velocity - target_velocity) < 1e10;
        assert(("Velocity matches the target velocity.", velocityMatchedTarget));

        bool yVelocityMatchedTarget = abs(ukf.x_(UKF_index::y) / delta_t - target_velocity) < 1e10;
        assert(("Y Velocity matches the target velocity.", yVelocityMatchedTarget));

        bool xVelocityIsZero = abs(ukf.x_(UKF_index::x) / delta_t) < 1e10;
        assert(("X Velocity is zero.", xVelocityIsZero));

        prev_x_position = ukf.x_(UKF_index::x);
        prev_y_position = ukf.x_(UKF_index::y);
    }

    bool averageVelocityMatchedTarget = abs((total_distance/total_t) - target_velocity) < 1e10;
    assert(("Average velocity across the whole distance matches the target velocity.", averageVelocityMatchedTarget));
    std::cout << "TestStraightLineConstantVelocity Completed\n";
};


// Checks that UKF predicted yaw rate is actually constant when the theta_acc model parameter is constant.
void TestConstantTurningRate() {
    UKF ukf;

    // Feed in first measurement
    MeasurementPackage meas_package;
    meas_package.sensor_type_ = MeasurementPackage::LASER;
    meas_package.raw_measurements_ = Eigen::VectorXd(2);
    meas_package.raw_measurements_ << -9.6649140872479418, 3.9175380731566825;
    meas_package.timestamp_ = 0;
    ukf.ProcessMeasurement(meas_package);

    // Set ukf theta_acc to target angular acceleration so we are turning at a constant rate
    double target_yaw_rate = M_PI/2.9107;
    ukf.x_(UKF_index::theta_acc) = target_yaw_rate;

    // Set ukf velocity to constant, non-zero rate
    ukf.x_(UKF_index::velocity) = -9.273;

    std::vector<double> delta_t_list = {1, 286, 77977};
    for (int i = 0; i < delta_t_list.size(); i++) {
        delta_t_list[i] = delta_t_list[i] / 1000000.0;
    }

    double prev_theta = ukf.x_(UKF_index::theta);

    for (auto& delta_t : delta_t_list) {
        ukf.Prediction(std::round(delta_t * 1000000));

        bool thetaMatchesTarget = abs(NormaliseAngle(ukf.x_(UKF_index::theta)) - NormaliseAngle(prev_theta + target_yaw_rate*delta_t)) < 1e10;
        assert(("Theta matches expected theta from constant yaw rate.", thetaMatchesTarget));

        prev_theta = ukf.x_(UKF_index::theta);
    }

    std::cout << "TestConstantTurningRate Completed\n";
};

void RunAllTests() {
    std::cout << "Running all tests\n";
    TestZeroDurationPredict();
    TestStraightLineConstantVelocity();
    TestConstantTurningRate();
    std::cout << "All tests Completed\n";
};

int main(int argc, char** argv) {
    RunAllTests();
}