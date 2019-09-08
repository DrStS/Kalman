#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  is_initialized_ = 0;

  // State dimension
  n_x_ = 5;
  n_z_ = 3;
  n_z_lidar_ = 2;
  // mean predicted measurement radar
  z_pred_ = VectorXd(n_z_);
  // measurement covariance matrix radar
  S_pred_ = MatrixXd(n_z_, n_z_);

  z_lidar_pred_ = VectorXd(n_z_lidar_);
  S_lidar_pred_ = MatrixXd(n_z_lidar_, n_z_lidar_);

  // Augmented state dimension
  n_aug_ = 7;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // create matrix for sigma points in measurement space
  Zsig_lidar_pred_ = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

  // create matrix for sigma points in measurement space
  Zsig_pred_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.90;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.50;

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
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  static unsigned int numCalls = 0;
  std::cout << "=========================" << std::endl;
  std::cout << "Current # function call UKF::ProcessMeasurement: " << numCalls << " time stamp: "<<meas_package.timestamp_<<std::endl;
  if (is_initialized_ < 2)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::LASER)
    {
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      std::cout << "Prediction for dT 0.0" << std::endl;
      P_ = MatrixXd::Identity(n_x_, n_x_) * (-0.99);
    }
    if (meas_package.sensor_type_ == MeasurementPackage::SensorType::RADAR)
    {
      x_(2) = (meas_package.raw_measurements_[2]) * cos(meas_package.raw_measurements_[1]);
      x_(3) = 0;
      x_(4) = 0;
    }  
    previous_time_stamp_ = meas_package.timestamp_;
    is_initialized_++;
  }

  static float dT = 0;
  if (previous_time_stamp_ != meas_package.timestamp_)
  {
    dT = (meas_package.timestamp_ - previous_time_stamp_) / 1000000.0;
    previous_time_stamp_ = meas_package.timestamp_;

    Prediction(dT);
    std::cout << "Prediction for dT " << dT << std::endl;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && meas_package.timestamp_>0)
  {
    std::cout << "UpdateRadar" << std::endl;
    UpdateRadar(meas_package);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::LASER && meas_package.timestamp_>0)
  {
    std::cout << "UpdateLidar" << std::endl;
    UpdateLidar(meas_package);
  }
  std::cout << "=========================" << std::endl;
  numCalls++;
//   if (numCalls > 14)
//    exit(EXIT_SUCCESS); //DEBUG
}

void UKF::Prediction(double delta_t)
{
  static unsigned int numCalls = 0;
  bool verbose = false;
  /**
   * Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  MatrixXd Xsig;
  GenerateSigmaPoints(&Xsig);
  MatrixXd Xsig_aug;
  AugmentedSigmaPoints(&Xsig_aug);
  SigmaPointPrediction(Xsig_aug, delta_t); // this writes memeber var Xsig_pred_
  PredictMeanAndCovariance();
  if (verbose)
  {
    std::cout.precision(5);
    std::cout << std::scientific;
    std::cout << "Current # function call UKF::Prediction: " << numCalls << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Xsig = " << Xsig << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Xsig_aug = " << std::endl
              << Xsig_aug << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Xsig_pred = " << std::endl
              << Xsig_pred_ << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Predicted state" << std::endl;
    std::cout << x_ << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Predicted covariance matrix" << std::endl;
    std::cout << P_ << std::endl;
    std::cout << "======================" << std::endl;
    //exit(EXIT_SUCCESS);//DEBUG
  }
  PredictRadarMeasurement();
  PredictLidarMeasurement();
  if (verbose)
  {
    std::cout << "======================" << std::endl;
    std::cout << "Predicted state radar" << std::endl;
    std::cout << z_pred_ << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Predicted sigma points radar" << std::endl;
    std::cout << Zsig_pred_ << std::endl;
    std::cout << "======================" << std::endl;
    std::cout << "Predicted covariance matrix radar" << std::endl;
    std::cout << S_pred_ << std::endl;
    std::cout << "======================" << std::endl;
  }

  numCalls++;
}

void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
   *  Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  bool verbose = false;
  VectorXd z = VectorXd(n_z_lidar_);
  z << meas_package.raw_measurements_[0], //px
      meas_package.raw_measurements_[1];  // py
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_lidar_pred_.col(i) - z_lidar_pred_;
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }
  // Kalman gain K;
  MatrixXd K = Tc * S_lidar_pred_.inverse();
  // residual
  VectorXd z_diff = z - z_lidar_pred_;
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_lidar_pred_ * K.transpose();
  double lidarNIS = z_diff.transpose() * S_lidar_pred_ * z_diff;
  // print result
  if (verbose)
    std::cout << "Time stamp: " << meas_package.timestamp_ << std::endl;
  std::cout << "Lidar NIS (Normalized Innovation Squared) value: " << lidarNIS << std::endl;
  if (verbose)
    std::cout << "Updated state x: " << std::endl
              << x_ << std::endl;
  //std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  // exit(EXIT_SUCCESS); //DEBUG
}

void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
   * Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  bool verbose = false;
  VectorXd z = VectorXd(n_z_);
  z << meas_package.raw_measurements_[0], // rho in m
      meas_package.raw_measurements_[1],  // phi in rad
      meas_package.raw_measurements_[2];  // rho_dot in m/s
  // create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  // calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_pred_.col(i) - z_pred_;
    // angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain K;
  MatrixXd K = Tc * S_pred_.inverse();
  // residual
  VectorXd z_diff = z - z_pred_;
  // angle normalization
  while (z_diff(1) > M_PI)
    z_diff(1) -= 2. * M_PI;
  while (z_diff(1) < -M_PI)
    z_diff(1) += 2. * M_PI;
  // update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_pred_ * K.transpose();

  double radarNIS = z_diff.transpose() * S_pred_ * z_diff;

  // print result
  if (verbose)
    std::cout << "Time stamp: " << meas_package.timestamp_ << std::endl;
  std::cout << "Radar NIS (Normalized Innovation Squared) value: " << radarNIS << std::endl;
  if (verbose)
    std::cout << "Updated state x: " << std::endl
              << x_ << std::endl;
}

//Private methods

void UKF::GenerateSigmaPoints(MatrixXd *Xsig_out)
{
  // define spreading parameter
  double lambda = 3 - n_x_;
  // create sigma point matrix
  MatrixXd Xsig = MatrixXd(n_x_, 2 * n_x_ + 1);
  // calculate square root of P
  MatrixXd A = P_.llt().matrixL();
  // set first column of sigma point matrix
  Xsig.col(0) = x_;
  // set remaining sigma points
  for (int i = 0; i < n_x_; ++i)
  {
    Xsig.col(i + 1) = x_ + sqrt(lambda + n_x_) * A.col(i);
    Xsig.col(i + 1 + n_x_) = x_ - sqrt(lambda + n_x_) * A.col(i);
  }
  // write result
  *Xsig_out = Xsig;
}

void UKF::AugmentedSigmaPoints(MatrixXd *Xsig_aug_out)
{
  // define spreading parameter
  lambda_ = 3 - n_aug_;
  // create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  // create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  // create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  // create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  // create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;
  // create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  // create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  for (int i = 0; i < n_aug_; ++i)
  {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }
  // write result
  *Xsig_aug_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd &Xsig_aug, double delta_t)
{
  // create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);
  // predict sigma points
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    // extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);
    // predicted state values
    double px_p, py_p;
    // avoid division by zero
    if (fabs(yawd) > 0.00001)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    }
    else
    {
      px_p = p_x + v * delta_t * cos(yaw);
      py_p = p_y + v * delta_t * sin(yaw);
    }
    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;
    // add noise
    px_p = px_p + 0.5 * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5 * nu_a * delta_t * delta_t * sin(yaw);
    v_p = v_p + nu_a * delta_t;
    yaw_p = yaw_p + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;
    // write predicted sigma point into right column
    Xsig_pred(0, i) = px_p;
    Xsig_pred(1, i) = py_p;
    Xsig_pred(2, i) = v_p;
    Xsig_pred(3, i) = yaw_p;
    Xsig_pred(4, i) = yawd_p;
  }
  // write result
  Xsig_pred_ = Xsig_pred;
}

void UKF::PredictMeanAndCovariance()
{
  // create vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  // set weights
  double weight_0 = lambda_ / (lambda_ + n_aug_);
  weights_(0) = weight_0;
  for (int i = 1; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 weights
    double weight = 0.5 / (n_aug_ + lambda_);
    weights_(i) = weight;
  }
  // predicted state mean
  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  // predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // iterate over sigma points
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    while (x_diff(3) > M_PI)
      x_diff(3) -= 2. * M_PI;
    while (x_diff(3) < -M_PI)
      x_diff(3) += 2. * M_PI;

    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }
}

void UKF::PredictRadarMeasurement(void)
{
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    // measurement model
    Zsig_pred_(0, i) = sqrt(p_x * p_x + p_y * p_y);                         // r
    Zsig_pred_(1, i) = atan2(p_y, p_x);                                     // phi
    Zsig_pred_(2, i) = (p_x * v1 + p_y * v2) / sqrt(p_x * p_x + p_y * p_y); // r_dot
  }

  // mean predicted measurement
  z_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_pred_ = z_pred_ + weights_(i) * Zsig_pred_.col(i);
  }

  // innovation covariance matrix S
  S_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_pred_.col(i) - z_pred_;

    // angle normalization
    while (z_diff(1) > M_PI)
      z_diff(1) -= 2. * M_PI;
    while (z_diff(1) < -M_PI)
      z_diff(1) += 2. * M_PI;

    S_pred_ = S_pred_ + weights_(i) * z_diff * z_diff.transpose();
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_, n_z_);
  R << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  S_pred_ = S_pred_ + R;
}

void UKF::PredictLidarMeasurement(void)
{
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // extract values for better readability
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    // measurement model
    Zsig_lidar_pred_(0, i) = p_x;
    Zsig_lidar_pred_(1, i) = p_y;
  }
  // mean predicted measurement
  z_lidar_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  {
    z_lidar_pred_ = z_lidar_pred_ + weights_(i) * Zsig_lidar_pred_.col(i);
  }
  // innovation covariance matrix S
  S_lidar_pred_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i)
  { // 2n+1 simga points
    // residual
    VectorXd z_diff = Zsig_lidar_pred_.col(i) - z_lidar_pred_;
    S_lidar_pred_ = S_lidar_pred_ + weights_(i) * z_diff * z_diff.transpose();
  }
  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
  R << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
  S_lidar_pred_ = S_lidar_pred_ + R;
}