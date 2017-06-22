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

  // Measurement dimension: radar measures r, phi, r_dot
  n_z_ = 3;
  
  // Sigma point spreading parameter
  lambda_ = 3 - n_aug_;

  // Weights of sigma points
  weights_ = VectorXd(2*n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i<2*n_aug_+1; i++) {  
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }

  // matrix with predicted sigma points in state space
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // matrix with sigma points in measurement space
  MatrixXd Zsig_ = MatrixXd(n_z_, 2 * n_aug_ + 1);

  // vector for mean predicted measurement
  VectorXd z_pred_ = VectorXd(n_z_);

  // matrix for predicted measurement covariance
  MatrixXd S_ = MatrixXd(n_z_, n_z_);

  // vector for incoming radar measurement
  VectorXd z_ = VectorXd(n_z_);

  // matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);

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
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x_(0) = meas_package.raw_measurements_[0];
      x_(1) = meas_package.raw_measurements_[1];
    }

    // done initializing
    previous_timestamp_ = meas_package.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
  *  Prediction & Update
  ****************************************************************************/
  float dt = (meas_package.timestamp_ - previous_timestamp_) / 1000000.0; //in secs
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    UpdateRadar(meas_package);

  } else {
    // Laser updates
    UpdateLidar(meas_package);
  }

  previous_timestamp_ = meas_package.timestamp_;

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
  VectorXd x_aug = VectorXd(7);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);

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
  MatrixXd x_mat = MatrixXd(n_aug, n_aug);
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


  /*****************************************************************************
  *  Predict Mean & Covariance
  ****************************************************************************/

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

  /*****************************************************************************
  *  Update Lidar Measurement
  ****************************************************************************/
  

  double NIS = z_diff.transpose() * S.inverse() * z_diff;
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

  /*****************************************************************************
  *  Predict Radar Sigma Points
  ****************************************************************************/


  /*****************************************************************************
  *  Update Radar
  ****************************************************************************/


  double NIS = z_diff.transpose() * S.inverse() * z_diff;
}
