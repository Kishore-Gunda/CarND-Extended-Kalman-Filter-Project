#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &Hj_in, MatrixXd &R_in,
                        MatrixXd &R_rad_in, MatrixXd &Q_in){
  cout<<"kalman init"<<endl;
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  H_radar_ = Hj_in;
  R_ = R_in;
  R_rad_ = R_rad_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state - updated from lecture
  */
 x_ = F_*x_;
 MatrixXd Ft = F_.transpose();
 P_ = F_*P_*Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations - updated form lecture
  */
  VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
  MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float p_x = x_[0];
  float p_y = x_[1];
  float v_x = x_[2];
  float v_y = x_[3];

  H_radar_ = tools.CalculateJacobian( x_ );
  VectorXd z_pred(3);
  float rho = sqrt( p_x*p_x + p_y*p_y );

  //return to avoid division by zero
  if( rho == 0.0)
    return;

  //create z vector in polar coordinates
  z_pred << rho, atan2( p_y, p_x ), ( p_x*v_x + p_y*v_y )/rho;

  // Update the state with ekf equations
  VectorXd y = z - z_pred;

  //Normalize angles 
  if( y[1] > M_PI )
    y[1] -= 2.0*M_PI;
  if( y[1] < -M_PI )
    y[1] += 2.0*M_PI;
  
  //Update final state - same as kalman filter
  MatrixXd H_radar_trans = H_radar_.transpose();
  MatrixXd P_Ht = P_*H_radar_trans;
  MatrixXd S = H_radar_* P_Ht + R_rad_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_Ht*Si;

  // Compute new state
  x_ = x_ + ( K*y );
  MatrixXd I_ = MatrixXd::Identity(x_.size(), x_.size());
  P_ = ( I_ - K*H_radar_)*P_;
}
