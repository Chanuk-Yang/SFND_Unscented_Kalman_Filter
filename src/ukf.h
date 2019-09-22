#ifndef UKF_H
#define UKF_H

#include <iostream>
#include "Eigen/Dense"
#include "measurement_package.h"

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

  // state covariance matrix
  Eigen::MatrixXd P_;

  // augmented state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate nu_a nu_yawdd]
  Eigen::VectorXd x_aug;

  // augmented state covariance matrix
  Eigen::MatrixXd P_aug;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

    // predicted sigma points matrix
  Eigen::MatrixXd Xsig_aug_;

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

  // Weight of sigma points at center
  double weight_0;

  // Weight of sigma points at distance
  double weight_i;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_aug_;

  // Sigma point spreading parameter
  double lambda_;

  // Lidar Measurement covariance matrix
  Eigen::MatrixXd R_lidar;

  // Lidar observation matrix
  Eigen::MatrixXd H_;

  // Lidar previous timestamp
  long lidar_timestamp_p;

  // Number of Radar measurement
  int n_z_radar;

  // Lidar Measurement covariance matrix
  Eigen::MatrixXd R_radar;

};

#endif  // UKF_H