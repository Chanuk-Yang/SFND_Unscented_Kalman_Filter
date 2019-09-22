#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

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

  // initial state vector
  x_aug = VectorXd(7);

  // initial covariance matrix
  P_aug = MatrixXd(7, 7);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
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
  
  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  x_ << 0,
        0,
        0,
        0,
        0;
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;
  R_ = MatrixXd(2,2);
  R_ << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  H_ = MatrixXd(2,5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;

  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  time_us_ = 0;

  Xsig_aug_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  weights_ = VectorXd(2*n_aug_ + 1);
  weight_0 = lambda_/(lambda_ + n_aug_);
  weight_i = 0.5/(lambda_ + n_aug_); 
  weights_ << weight_0,
              weight_i, weight_i, weight_i, weight_i, weight_i, weight_i, weight_i, 
              weight_i, weight_i, weight_i, weight_i, weight_i, weight_i, weight_i;
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  if (meas_package.sensor_type_ == meas_package.LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == meas_package.RADAR)
  {
    UpdateRadar(meas_package);
  }
  else
  {
    std::cout << "MEASUREMEENT ERROR, CHECK IF THE DATA IS LASER OR RADAR!" <<std::endl;
  }
  
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location. 
   * Modify the state vector, x_. Predict sigma points, the state, 
   * and the state covariance matrix.
   */
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd SRM = P_aug.llt().matrixL();

  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; ++i) {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * SRM.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * SRM.col(i);
  }

  for (int i = 0; i< 2*n_aug_+1; ++i) 
  {
    double p_x = Xsig_aug_(0,i);
    double p_y = Xsig_aug_(1,i);
    double v = Xsig_aug_(2,i);
    double yaw = Xsig_aug_(3,i);
    double yawd = Xsig_aug_(4,i);
    double nu_a = Xsig_aug_(5,i);
    double nu_yawdd = Xsig_aug_(6,i);

    double px_pred, py_pred;

    if (fabs(yawd) > 0.0001) 
    {
        px_pred = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_pred = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    } 
    else 
    {
        px_pred = p_x + v*delta_t*cos(yaw);
        py_pred = p_y + v*delta_t*sin(yaw);
    }

    double v_pred = v;
    double yaw_pred = yaw + yawd*delta_t;
    double yawd_pred = yawd;

    px_pred = px_pred + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_pred = py_pred + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_pred = v_pred + nu_a*delta_t;

    yaw_pred = yaw_pred + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_pred = yawd_pred + nu_yawdd*delta_t;

    Xsig_pred_(0,i) = px_pred;
    Xsig_pred_(1,i) = py_pred;
    Xsig_pred_(2,i) = v_pred;
    Xsig_pred_(3,i) = yaw_pred;
    Xsig_pred_(4,i) = yawd_pred;
  }

  VectorXd x_sum = VectorXd(n_x_);
  x_sum.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  
    x_sum = x_sum + weights_(i) * Xsig_pred_.col(i);
  }

  MatrixXd P_sum = MatrixXd(n_x_,n_x_);
  P_sum.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  

    VectorXd x_diff = Xsig_pred_.col(i) - x_sum;

    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

    P_sum = P_sum + weights_(i) * x_diff * x_diff.transpose();
  }

  x_ = x_sum;
  P_ = P_sum;
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */

  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);

  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
}