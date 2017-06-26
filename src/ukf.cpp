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
  P_ = MatrixXd::Identity(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // fastest measured linear accel for street legal car is 
  // 0-60 mph in 2.2 second (12 m/s^2)
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

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
  TO_DID:
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  // State dimension
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;
  
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  double w = 0.5/(n_aug_+lambda_);
  for (int i=1; i<2*n_aug_+1; i++) {  
    weights_(i) = w;
  }

  // matrix with predicted sigma points in state space
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TO_DID:
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */

  /*****************************************************************************
  *  Initialization
  ****************************************************************************/
  long previous_timestamp;

  if (!is_initialized_) {
    // first measurement
    cout << "UKF: " << endl;
    x_.fill(0.0);  // values are important to RMSE

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];
      float rho_dot = meas_package.raw_measurements_[2];

      x_(0) = rho * cos(phi);
      x_(1) = rho * sin(phi);
      cout << "x_: " << x_ << endl;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
      cout << "x_: " << x_ << endl;
    }

    // done initializing
    previous_timestamp = meas_package.timestamp_;
    cout << "previous_timestamp: " << previous_timestamp << endl;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
  *  Prediction & Update
  ****************************************************************************/
  float delta_t = (meas_package.timestamp_ - previous_timestamp) / 1000000.0; //in secs
  cout << "delta_t: " << delta_t << endl;
  Prediction(delta_t);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);

  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  previous_timestamp = meas_package.timestamp_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODOING:
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  /*****************************************************************************
  *  Create Augmented Sigma Points
  ****************************************************************************/

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_ + 1);

  //create augmented mean state
  x_aug.head(n_x_) = x_;
  x_aug(n_aug_-2) = 0;
  x_aug(n_aug_-1) = 0;

  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(n_aug_-2, n_aug_-2) = pow(std_a_, 2);
  P_aug(n_aug_-1, n_aug_-1) = pow(std_yawdd_ ,2);

  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd x_mat = MatrixXd(n_aug_, n_aug_);
  x_mat << x_aug, x_aug, x_aug, x_aug, x_aug, x_aug, x_aug;
  //std::cout << x_mat << std::endl;
  
  //calculate sigma points ...
  //set sigma points as columns of matrix Xsig
  Xsig_aug.col(0) = x_aug;
  
  Xsig_aug.block(0,1, n_aug_,n_aug_) = x_mat + (sqrt(lambda_+n_aug_) * A_aug);
  //std::cout << Xsig_aug.block(0,1, n_x,n_x) << std::endl;
  
  Xsig_aug.block(0,n_aug_+1, n_aug_,n_aug_) = x_mat - (sqrt(lambda_+n_aug_) * A_aug);

  /*****************************************************************************
  *  Predict Sigma Points
  ****************************************************************************/

  //placeholder for sigma points predictions
  VectorXd noise = VectorXd(n_x_);
  VectorXd preds = VectorXd(n_x_);

  //loop through sigma points columns
  for (int i=0; i<(2*n_aug_+1); i++) {

    // assign values
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    noise << 0.5 * pow(delta_t,2) * cos(yaw) * nu_a,
        0.5 * pow(delta_t,2) * sin(yaw) * nu_a,
        delta_t * nu_a,
        0.5 * pow(delta_t,2) * nu_yawdd,
        delta_t * nu_yawdd;
    
    //avoid division by zero
    if (fabs(yawd) < 0.0001) {
      
      //std::cout << "Div by zero:" << std::endl << Xsig_aug.col(i) << std::endl;
      preds << v * cos(yaw) * delta_t,
          v * sin(yaw) * delta_t,
          0, yawd*delta_t, 0;

      preds = Xsig_aug.block(0,i, 5,1) + preds + noise;

    } else {

      preds << (v/yawd) * ( sin(yaw + yawd*delta_t) - sin(yaw) ),
          (v/yawd) * ( -cos(yaw + yawd*delta_t) + cos(yaw) ),
          0, yawd*delta_t, 0;

      preds = Xsig_aug.block(0,i, 5,1) + preds + noise;
    }
    
    //write predicted sigma points into right column
    Xsig_pred_.col(i) = preds;
  }

  /*****************************************************************************
  *  Predict Mean & Covariance
  ****************************************************************************/

  //predict state mean
  x_.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i=0; i<2*n_aug_ + 1; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;  //state diff

    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.0*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.0*M_PI;

    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODOING:
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  /*****************************************************************************
  *  Update Lidar Measurement
  ****************************************************************************/
  
  //create vector for incoming lidar measurement
  VectorXd z = meas_package.raw_measurements_;

  //measurement matrix - laser
  MatrixXd H = MatrixXd(2, n_x_);
  H << 1, 0, 0, 0, 0, 
      0, 1, 0, 0, 0;

  //measurement covariance matrix - laser
  MatrixXd R = MatrixXd(2, 2);
  R << std_laspx_*std_laspx_, 0, 
      0, std_laspy_*std_laspy_;

  VectorXd y = z - H * x_;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  MatrixXd I = MatrixXd::Identity(n_x_, n_x_);
  P_ = (I - (K * H)) * P_;

  //calculate Normalized Innovation Square (NIS)
  //double NIS = y.transpose() * Si * y;
  cout << "Lidar update, x_: " << x_ << endl;
  cout << "Lidar update, P_: " << P_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODOING:
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  /*****************************************************************************
  *  Predict Radar Sigma Points
  ****************************************************************************/

  // Measurement dimension: radar measures r, phi, r_dot
  int n_z = 3;

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);

  //transform sigma points into measurement space
  Zsig.fill(0.0);
  for (int i=0; i<2*n_aug_ + 1; i++) {
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);
    double yawd = Xsig_pred_(4,i);

    float zero_check = 0.0001;
    if (fabs(p_x)<0.0001) {
      if (fabs(p_y)<0.0001) {
        p_y = zero_check;
      }
      p_x = zero_check;
    }

    //measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);  //r
    Zsig(1,i) = atan(p_y/p_x);  //phi
    Zsig(2,i) = (p_x*cos(yaw)*v + p_y*sin(yaw)*v) / Zsig(0,i);  //r_dot
  }
  
  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i=0; i<2*n_aug_ + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;  //residual

    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  MatrixXd R = MatrixXd(n_z,n_z);
  R(0,0) = std_radr_*std_radr_;
  R(1,1) = std_radphi_*std_radphi_;
  R(2,2) = std_radrd_*std_radrd_;

  S = S + R;

  /*****************************************************************************
  *  Update Radar
  ****************************************************************************/

  //create vector for incoming radar measurement
  VectorXd z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i=0; i<2*n_aug_+1; i++) {    
    VectorXd z_diff = Zsig.col(i) - z_pred;  //residual
    
    //angle normalization
    while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
    while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;

    VectorXd x_diff = Xsig_pred_.col(i) - x_;  // state difference
    
    //angle normalization
    while (x_diff(3) > M_PI) x_diff(3) -= 2.0*M_PI;
    while (x_diff(3) < -M_PI) x_diff(3) += 2.0*M_PI;

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  VectorXd z_diff = z - z_pred;  //residual

  //angle normalization
  while (z_diff(1) > M_PI) z_diff(1) -= 2.0*M_PI;
  while (z_diff(1) < -M_PI) z_diff(1) += 2.0*M_PI;  

  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //calculate Normalized Innovation Square (NIS)
  //double NIS = z_diff.transpose() * S.inverse() * z_diff;
  cout << "Radar update, x_: " << x_ << endl;
  cout << "Radar update, P_: " << P_ << endl;
}
