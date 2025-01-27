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
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 3;
  
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

  // state initialize
  is_initialized_ = false;
  n_x_ = 5;
  n_aug_ = 7;
  time_us_ = 0;

  x_ << 0,
        0,
        0,
        0,
        0;
  P_ << std_laspx_*std_laspx_, 0, 0, 0, 0,
        0, std_laspy_*std_laspy_, 0, 0, 0,
        0, 0, std_radrd_*std_radrd_, 0, 0,
        0, 0, 0, std_radphi_*std_radrd_, 0,
        0, 0, 0, 0, std_radphi_*std_radphi_;

  // Lidar Configuration
  R_lidar = MatrixXd(2,2);
  R_lidar << std_laspx_*std_laspx_, 0,
        0, std_laspy_*std_laspy_;
  H_ = MatrixXd(2,5);
  H_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0;
  lidar_timestamp_p = 0;

  // Radar configuration
  n_z_radar = 3;
  R_radar = MatrixXd(n_z_radar,n_z_radar);
  R_radar <<  std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

  // UKF configuration
  lambda_ = 3 - n_aug_;
  Xsig_aug_ = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

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
  if (!is_initialized_)
  {
    if (meas_package.sensor_type_ == meas_package.LASER)
    {
      x_(0) = meas_package.raw_measurements_(0);
      x_(1) = meas_package.raw_measurements_(1);
      // x_(2) = 0;
      // x_(3) = 0.0;
      // x_(4) = 0.0;
      is_initialized_ = true;
      // std::cout << "init lidar" << std::endl;

    }
    else if (meas_package.sensor_type_ == meas_package.RADAR)
    {
      double r = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rd = meas_package.raw_measurements_(2);
      x_(0) = r * cos(phi);
      x_(1) = r * sin(phi);
      x_(2) = rd;
      x_(3) = phi;
      // x_(4) = 0.0;
      is_initialized_ = true;
      std::cout << "init radar" << std::endl;
    }
    time_us_ = meas_package.timestamp_;
  }
  else
  {
    // std::cout << meas_package.timestamp_ <<" and "<< time_us_ << std::endl;
    double dt = static_cast<double>(meas_package.timestamp_ - time_us_) / 1e6;
    // std::cout << meas_package.sensor_type_ << " dt is "<<dt << std::endl;
    if (meas_package.sensor_type_ == meas_package.LASER)
    {
      Prediction(dt);
      UpdateLidar(meas_package);
      time_us_ = meas_package.timestamp_;
    }
    else if (meas_package.sensor_type_ == meas_package.RADAR)
    {
      Prediction(dt);
      UpdateRadar(meas_package);
      time_us_ = meas_package.timestamp_;
    }
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

  Xsig_pred_.setZero();
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

    if (fabs(yawd) > 0.00000001) 
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
  
  lidar_timestamp_p = meas_package.timestamp_;
  VectorXd z = VectorXd(2);
  z << meas_package.raw_measurements_;
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_lidar;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  x_ = x_ + (K * y);

  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - K * H_) * P_;
  
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  // create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar, 2 * n_aug_ + 1);

  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar);
  

  // measurement model predict in each sigma point
  
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  { 
    
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                   
    Zsig(1,i) = atan2(p_y,p_x);                       
    Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y); 
  }

  // sum with weighed sigma point
  z_pred.setZero();
  for (int i=0; i < 2*n_aug_+1; ++i) 
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }
  // innovate covariance matrix and calculate kalman gain
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar);
  MatrixXd S = MatrixXd(n_z_radar,n_z_radar);
  Tc.setZero();
  S.setZero();
  for (int i = 0; i < 2 * n_aug_ + 1; ++i) 
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;


    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // // angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    S  += weights_(i) * z_diff * z_diff.transpose();
    Tc += weights_(i) * x_diff * z_diff.transpose();
  }
  S += R_radar;
  MatrixXd K = Tc * S.inverse();

  // calculate z difference
  VectorXd z = VectorXd(n_z_radar);
  z << meas_package.raw_measurements_;
  VectorXd z_dev = z - z_pred;

  // angle normalization
  while (z_dev(1)> M_PI) z_dev(1)-=2.*M_PI;
  while (z_dev(1)<-M_PI) z_dev(1)+=2.*M_PI;

  // update state mean and covariance matrix
  x_ += K * z_dev;
  P_ -= K*S*K.transpose();
}