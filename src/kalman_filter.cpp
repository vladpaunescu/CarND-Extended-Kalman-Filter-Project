#include "kalman_filter.h"
#include <cmath>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
		MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict() {
	// kalman filter state prediction

	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
	/**
  TODO:
	 * update the state by using Kalman Filter equations
	 */

	// did refactor common state updates from
	// normal Kalman Filter and Extended Kalman Filter
	// into one internal method

	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;

	// internal update
	UpdateInternal(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
	/**
  TODO:
	 * update the state by using Extended Kalman Filter equations
	 *
	 */

	VectorXd z_pred = h(x_);
	VectorXd y = z - z_pred;

	// phi normalization
	y[1] = y[1] < -M_PI ? M_PI + fmod(y[1] + M_PI, 2 * M_PI) : fmod(y[1] + M_PI, 2 * M_PI) - M_PI;

	// internal update
	UpdateInternal(y);
}

// private stuff

VectorXd KalmanFilter::h(const VectorXd& predicted_state) {
	double px = predicted_state[0];
	double py = predicted_state[1];
	double vx = predicted_state[2];
	double vy = predicted_state[3];

	VectorXd projected_state = VectorXd(3);
	projected_state << 0, 0, 0;

	if(px == 0 && py == 0) {
		return projected_state;
	}

	double rho = sqrt(px * px + py * py);
	double phi = atan2(py, px);
	double rho_dot = ((px * vx) + (py * vy)) / rho;

	projected_state << rho, phi, rho_dot;

	return projected_state;
}

// common update step given the error
void KalmanFilter::UpdateInternal(const VectorXd& y) {
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
