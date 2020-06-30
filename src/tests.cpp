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
    std::cout << "Starting TestZeroDurationPredict...\n";
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
//    std::cout << "x_diff=\n" << x_diff << "\n";
//    std::cout << "P_diff=\n" << P_diff << "\n";

    bool x_mean_changed = x_diff.isZero(1e-7);
    assert(("Mean state unchanged between initialization and zero-duration predict step.", x_mean_changed));

    bool P_mean_changed = P_diff.isZero(1e-7);
    assert(("Mean covariance matrix unchanged between initialization and zero-duration predict step.", P_mean_changed));

    // Store UKF state
    x_prev = ukf.x_;
    P_prev = ukf.P_;

    // Predict state after 0 time
    ukf.Prediction(0);

    x_diff = x_prev - ukf.x_;
    P_diff = P_prev - ukf.P_;

    x_mean_changed = x_diff.isZero(1e-10);
    assert(("Mean state unchanged between first and second zero-duration predict step.", x_mean_changed));

    P_mean_changed = P_diff.isZero(1e-10);
    assert(("Mean covariance matrix unchanged between first and second zero-duration predict step.", P_mean_changed));
    std::cout << "TestZeroDurationPredict Completed\n";
};


// Checks that UKF predicted velocity is constant if we're travelling in a straight line.
void TestStraightLineConstantVelocity() {
    std::cout << "Starting TestStraightLineConstantVelocity...\n";
    UKF ukf;

    // Set ukf velocity to 5.81 m/s
    double target_velocity = 5.81;

    // Set up initial state
    ukf.x_(UKF_index::x) = 7.21;
    ukf.x_(UKF_index::y) = 2.86;
    ukf.x_(UKF_index::velocity) = target_velocity;
    ukf.x_(UKF_index::theta) = M_PI/2.0;
    ukf.x_(UKF_index::theta_acc) = 0.0;
//    ukf.P_ *= 1e-12;
    ukf.P_.fill(0.0);
    ukf.std_a_ = 0.0;
    ukf.std_yawdd_ = 0.0;
    ukf.is_initialized_ = true;

    // Set ukf theta_acc to zero so we are travelling in a straight line
    ukf.x_(UKF_index::theta_acc) = 0;

    std::vector<double> delta_t_list = {1, 168, 50000};
    int total_t = 0;
    for (auto& delta_t : delta_t_list) {
        total_t += delta_t;
    }
    double total_distance = 0;

    double prev_x_position = ukf.x_(UKF_index::x);
    double prev_y_position = ukf.x_(UKF_index::y);

    for (int i = 0; i < delta_t_list.size(); i++) {
        int delta_t = delta_t_list[i];
        ukf.Prediction(delta_t_list[i]);
        double x_diff = ukf.x_(UKF_index::x) - prev_x_position;
        double y_diff = ukf.x_(UKF_index::y) - prev_y_position;
        double distance = std::sqrt(x_diff * x_diff + y_diff * y_diff);
        total_distance += distance;

        double current_velocity = distance / delta_t;
        bool velocityMatchedTarget = abs(current_velocity - target_velocity) < 1e-6;
        assert(("Velocity matches the target velocity.", velocityMatchedTarget));

        double y_velocity = y_diff / delta_t;
        bool yVelocityMatchedTarget = abs(y_velocity - target_velocity) < 1e-6;
        assert(("Y Velocity matches the target velocity.", yVelocityMatchedTarget));

        double x_velocity = x_diff / delta_t;
        bool xVelocityIsZero = abs(x_velocity) < 1e-6;
        assert(("X Velocity is zero.", xVelocityIsZero));

        prev_x_position = ukf.x_(UKF_index::x);
        prev_y_position = ukf.x_(UKF_index::y);
    }

    double average_velocity = total_distance / total_t;
    bool averageVelocityMatchedTarget = abs(average_velocity - target_velocity) < 1e-7;
    assert(("Average velocity across the whole distance matches the target velocity.", averageVelocityMatchedTarget));
    std::cout << "TestStraightLineConstantVelocity Completed\n";
};


// Checks that UKF predicted yaw rate is actually constant when the theta_acc model parameter is constant.
void TestConstantTurningRate() {
    std::cout << "Starting TestConstantTurningRate...\n";
    UKF ukf;

    // Set ukf theta_acc to target angular acceleration so we are turning at a constant rate
    double target_yaw_rate = M_PI/2.9107;

    // Set up initial state
    ukf.x_(UKF_index::x) = 57.1;
    ukf.x_(UKF_index::y) = -6.29;
    ukf.x_(UKF_index::velocity) = 27;
    ukf.x_(UKF_index::theta) = M_PI * 0.64;
    ukf.x_(UKF_index::theta_acc) = target_yaw_rate;
//    ukf.P_ *= 1e-12;
    ukf.P_.fill(0.0);
    ukf.std_a_ = 0.0;
    ukf.std_yawdd_ = 0.0;
    ukf.is_initialized_ = true;

    ukf.x_(UKF_index::theta_acc) = target_yaw_rate;

    // Set ukf velocity to constant, non-zero rate
    ukf.x_(UKF_index::velocity) = -9.273;

    std::vector<double> delta_t_list = {1, 286, 77977};
    double prev_theta = ukf.x_(UKF_index::theta);

    for (int i = 0; i < delta_t_list.size(); i++) {
        double delta_t = delta_t_list[i];
        ukf.Prediction(delta_t_list[i]);

        double normalisedTheta = NormaliseAngle(ukf.x_(UKF_index::theta));
        double normalisedThetaTarget = NormaliseAngle(prev_theta + target_yaw_rate * delta_t);
        bool thetaMatchesTarget = abs(normalisedTheta - normalisedThetaTarget) < 1e-2;
        assert(("Theta matches expected theta from constant yaw rate.", thetaMatchesTarget));

        prev_theta = ukf.x_(UKF_index::theta);
    }

    std::cout << "TestConstantTurningRate Completed\n";
};

void RunAllTests() {
    std::cout << "##################### Running all tests\n";
    TestZeroDurationPredict();
    TestStraightLineConstantVelocity();
    TestConstantTurningRate();
    std::cout << "##################### All tests Completed\n";
};

int main(int argc, char** argv) {
    RunAllTests();
}