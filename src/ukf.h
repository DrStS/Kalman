#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"

class UKF
{
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

private:
  /**
   * Generate Sigma Point
   * @param[out] Xsig_out 
   */
  void GenerateSigmaPoints(Eigen::MatrixXd *Xsig_out);
  /**
   * Augmented Sigma Points
   * @param[out] Xsig_aug_out 
   */
  void AugmentedSigmaPoints(Eigen::MatrixXd *Xsig_aug_out);
  /**
   * Predict Sigma Points
   * @param[in] Xsig_aug
   * @param[in] delta_t
   */
  void SigmaPointPrediction(Eigen::MatrixXd &Xsig_aug, double delta_t);
  /**
   * Predict Mean and Covariance 
   */
  void PredictMeanAndCovariance(void);
  /**
   * Predict Radar Measurement
   */
  void PredictRadarMeasurement(void);
  /**
   * Predict Lidar Measurement 
   */
  void PredictLidarMeasurement(void);

public:
  // initially set to false, set to true in first call of ProcessMeasurement
  unsigned int is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state lidar vector:
  Eigen::VectorXd z_lidar_pred_;

  // state radar vector:
  Eigen::VectorXd z_pred_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // state covariance matrix in lidar measurement space
  Eigen::MatrixXd S_lidar_pred_;

  // state covariance matrix in radar measurement space
  Eigen::MatrixXd S_pred_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // predicted sigma points matrix in lidar measurement space
  Eigen::MatrixXd Zsig_lidar_pred_;

  // predicted sigma points matrix in radar measurement space
  Eigen::MatrixXd Zsig_pred_;

  // previous time stamp
  long previous_time_stamp_;

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
  double std_radrd_;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // State dimension radar
  int n_z_;

  // State dimension radar
  int n_z_lidar_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;
};

#endif // UKF_H