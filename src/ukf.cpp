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
  x_.fill(0.0);

  // initial covariance matrix
  P_ = MatrixXd::Identity(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.5; // 5 - 2

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.35; // 0.5 - 0.2

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

  is_initialized_ = false;
  time_us_ = 0.0;

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
  Xsig_pred_.fill(0.0);

  // Augmented Sigma Points Matrix
  Xsig_aug_ = MatrixXd(n_aug_, n_aug_sigma_);
  Xsig_aug_.fill(0.0);

  // Sigma Point Weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_.fill(0.0);

  // Set Sigma Point Weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;

  for (int i = 1; i < n_aug_sigma_; i++)
  {
    double weight = 0.5 /(n_aug_ + lambda_);
    weights_(i) = weight;
  }

  NIS_radar_ = 0.0;
  NIS_lidar_ = 0.0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  if (!is_initialized_) {
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = meas_package.raw_measurements_[0];
      double phi = meas_package.raw_measurements_[1];
      double rho_dot = meas_package.raw_measurements_[2];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      // x_ << px, py, 0, 0, 0;

      x_[0] = px;
      x_[1] = py;
      x_[2] = rho_dot * py;
      x_[3] = rho_dot * px;

    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      double px = meas_package.raw_measurements_[0];
      double py = meas_package.raw_measurements_[1];
      // x_ << px, py, 0, 0, 0;

      x_[0] = px;
      x_[1] = py;
    }

    if (fabs(x_[0]) < 0.001 && fabs(x_[1]) < 0.001 ) {
      x_[0] = 0.001;
      x_[1] = 0.001;
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

  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_ == true) {
    n_z_ = 3;
    UKF::UpdateRadar(meas_package);
    cout << "Update Radar" << endl;
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_ == true) {
    n_z_ = 2;
    UKF::UpdateLidar(meas_package);
    cout << "Update Lidar" << endl;
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
    Prediction:
      1. Generate Sigma Points
      2. Predict Sigma Points
      3. Predict Mean and Covariance
  */

  //* 1. Generate Sigma points
  UKF::GenerateSigmaPoints();

  //* 2. Predict Sigma Points
  UKF::SigmaPointPrediction(delta_t);

  //* 3. Predict Mean and Covariance
  UKF::PredictMeanAndCovariance();
}

void UKF::GenerateSigmaPoints() {
  // Create Augmented Mean Vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);

  // Create Augmented State Covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);

  // Xsig_aug is already delcared: Xsig_aug_

  // Create Augmented Mean State
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // Create Augmented Covariance Matrix
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_aug_-2, n_aug_-2) = std_a_ * std_a_;
  P_aug(n_aug_-1, n_aug_-1) = std_yawdd_ * std_yawdd_;

  // Create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  L.fill(0.0);

  // Create Augmented Sigma Points
  Xsig_aug_.col(0) = x_aug;
  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug_.col(i + 1)           = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug_.col(i + 1 + n_aug_)  = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // cout << "Xsig_aug_: " << Xsig_aug_ << endl;
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

void UKF::PredictMeanAndCovariance() {
  // !weights set in initialization: weights_

  // Predicted State Mean
  x_.fill(0.0);
  for (int i =0; i < n_aug_sigma_; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }

  // Predicted State Covariance
  P_.fill(0.0);
  for (int i =0; i < n_aug_sigma_; i++) {
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_ = P_ * weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  z_ = VectorXd(n_z_); // z_pred in assignment
  z_.fill(0.0);
  z_ << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1],
        meas_package.raw_measurements_[2];

  UKF::PredictRadarMeasurement();
  UKF::UpdateState(meas_package);

  NIS_radar_ = (z_ - z_pred_).transpose() * S_.inverse() * (z_ - z_pred_);
  cout << "NIS_radar_: " << NIS_radar_ << endl;
}

void UKF::PredictRadarMeasurement() {
  ///* Matrix for Sigma Points into Measurement Space
  Zsig_ = MatrixXd(n_z_, n_aug_sigma_);
  Zsig_.fill(0.0);

  // Transform sigma points into measurement space
  for (int i = 0; i < n_aug_sigma_; i++) {
    double p_x, p_y, v, yaw;

     p_x = Xsig_pred_(0, i);
     p_y = Xsig_pred_(1, i);
     v = Xsig_pred_(2, i);
     yaw = Xsig_pred_(3, i);

     double v1, v2;
     v1 = cos(yaw) * v;
     v2 = sin(yaw) * v;

    //  Measurement Model
    if (p_x == 0 && p_y == 0) {
      Zsig_(0, i) = 0;
      Zsig_(1, i) = 0;
      Zsig_(2, i) = 0;
    } else {
      Zsig_(0, i) = sqrt(p_x*p_x + p_y*p_y);
      Zsig_(1, i) = atan2(p_y, p_x);
      Zsig_(2, i) = (p_x * v1 + p_y * v2)/ sqrt(p_x*p_x + p_y*p_y);
    }
  }

  // Predicted Measurement Mean: z_pred
  z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < n_aug_sigma_; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  // Predicted Measurement Covariance
  S_ = MatrixXd(n_z_, n_z_);
  S_.fill(0.0);
  for (int i = 0; i < n_aug_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  // Add Measurement Noise Covariance Matrix
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R.fill(0.0);
  R <<  std_radr_*std_radr_, 0, 0,
        0, std_radphi_*std_radphi_, 0,
        0, 0, std_radrd_*std_radrd_;

  S_ = S_ + R;

  // cout << "z_pred_: " << endl << z_pred_ << endl;
  // cout << "S_: " << endl << S_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  z_ = VectorXd(n_z_); // z_pred in assignment
  z_.fill(0.0);
  z_ << meas_package.raw_measurements_[0],
        meas_package.raw_measurements_[1];

  UKF::PredictLidarMeasurement();
  UKF::UpdateState(meas_package);

  NIS_lidar_ = (z_ - z_pred_).transpose() * S_.inverse() * (z_ - z_pred_);
  cout << "NIS_lidar_: " << NIS_lidar_ << endl;
}

void UKF::PredictLidarMeasurement() {
  ///* Matrix for Sigma Points into Measurement Space
  Zsig_ = MatrixXd(n_z_, n_aug_sigma_);
  Zsig_.fill(0.0);

  // Transform sigma points into measurement space
  for (int i = 0; i < n_aug_sigma_; i++) {
    double p_x, p_y, v, yaw;

     p_x = Xsig_pred_(0, i);
     p_y = Xsig_pred_(1, i);

    //  Measurement Model
    Zsig_(0, i) = p_x;
    Zsig_(1, i) = p_y;
  }

  // Predicted Measurement Mean: z_pred
  z_pred_ = VectorXd(n_z_);
  z_pred_.fill(0.0);
  for (int i = 0; i < n_aug_sigma_; i++) {
    z_pred_ = z_pred_ + weights_(i) * Zsig_.col(i);
  }

  // Predicted Measurement Covariance
  S_ = MatrixXd(n_z_, n_z_);
  S_.fill(0.0);
  for (int i = 0; i < n_aug_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;

    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1)<-M_PI) z_diff(1) += 2. * M_PI;

    S_ = S_ + weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z_, n_z_);
  R.fill(0.0);
  R <<  std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;

  S_ = S_ + R;

  // cout << "z_pred_: " << endl << z_pred_ << endl;
  // cout << "S_: " << endl << S_ << endl;
}


void UKF::UpdateState(MeasurementPackage meas_package) {
  // Matrix for Cross Correlation
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  Tc.fill(0.0);

  for (int i=0; i< n_aug_sigma_; i++) {
    // residual
    VectorXd z_diff = Zsig_.col(i) - z_pred_;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // angle normalization: z_diff & x_diff
      while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
      while (z_diff(1) <-M_PI) z_diff(1) += 2. * M_PI;

      while (x_diff(3) > M_PI) x_diff(3) -= 2. * M_PI;
      while (x_diff(3) <-M_PI) x_diff(3) += 2. * M_PI;
    }

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman Gain: K
  MatrixXd K = Tc * S_.inverse();

  // residual, again
  VectorXd z_diff = z_ - z_pred_;

  // angle normalization: z_diff
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    while (z_diff(1) > M_PI) z_diff(1) -= 2. * M_PI;
    while (z_diff(1) <-M_PI) z_diff(1) += 2. * M_PI;
  }

  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_ * K.transpose();
}
