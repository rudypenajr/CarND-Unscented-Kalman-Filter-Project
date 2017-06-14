#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // State Dimension
  n_x_ = 5;

  // Augmented State Dimension
  n_aug_ = 7;

  // # of Sigma Points
  n_aug_sigma_ = 2 * n_aug_ + 1;

  // Sigma Point Spreading Parameter
  lambda_ = 3 - n_aug_;

  // Predicted Sigma Points Matrix
  Xsig_pred_ = MatrixXd(n_x_, n_aug_sigma_);

  // Augmented Sigma Points Matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_aug_sigma_);

  // Sigma Point Weights
  weights_ = VectorXd(2 * n_aug_ + 1);

  // Set Sigma Point Weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < n_aug_sigma_; i++)
  {
    double weight = 0.5 /(n_aug_ + lambda_);
    weights_(i) = weight;
  }

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      x_ << px, py, 0, 0, 0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      x_ << px, py, 0, 0, 0;
    }

    time_us_ = meas_package.timestamp_;
    is_initialized_ = true;

    cout << "Initialition Complete." << endl;
    cout << "x_: " << "\n" << x_ << endl;
    cout << "P_: " << "\n" << P_ << endl;
    return;
  }

  double delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  cout << "time_us_: " << time_us_ << endl;

  UKF::Prediction(delta_t);
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /**
    Prediction:
      1. Generate Sigma Points
      2. Predict Sigma Points
      3. Predict Mean and covariance
  */

  //* 1. Generate Sigma points
  UKF::GenerateSigmaPoints();

  //* 2. Predict Sigma Points
  UKF::SigmaPointPrediction(delta_t);
}

void UKF::GenerateSigmaPoints() {
  // Create Augmented Mean Vector
  VectorXd x_aug = VectorXd(7);

  // Create Augmented State Covariance
  MatrixXd P_aug = MatrixXd(7, 7);

  // Xsig_aug is already delcared: Xsig_aug_

  // Create Augmented Mean State
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create Augmented Covariance Matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_*std_a_;
  P_aug(6, 6) = std_yawdd_*std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  // Create Augmented Sigma Points
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_sigma_; i++) {
    Xsig_aug_.col(i + 1)           = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  cout << "Xsig_aug_: " << Xsig_aug_ << endl;
}

void UKF::SigmaPointPrediction(double delta_t) {
   for (int i = 0; i < n_aug_sigma_; i++) {
     double p_x, p_y, v, yaw, yawd, nu_a, nu_yawdd;

      p_x = Xsig_aug_(0, i);
      p_y = Xsig_aug_(1, i);
      v = Xsig_aug_(2, i);
      yaw = Xsig_aug_(3, i);
      yawd = Xsig_aug_(4, i);
      nu_a = Xsig_aug_(5, i);
      nu_yawdd = Xsig_aug_(6, i);

    double px_p, py_p; // predicted state values

    if (fabs(yawd) > 0.001) {
      px_p = p_x + v/yawd * ( sin(yaw + yawd*delta_t) - sin(yaw) );
      py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw + yawd*delta_t) );
    } else {
      px_p = p_x + v*delta_t*cos(yaw);
      py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
}
