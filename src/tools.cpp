#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  
  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() == 0 || estimations.size() != ground_truth.size())
  {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }
  
  //accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i){
    VectorXd est = estimations[i];
    VectorXd act = ground_truth[i];
    VectorXd diff = act - est;
    VectorXd mul = diff.array() * diff.array();
    rmse += mul;
  }
  
  //calculate the mean
  rmse /= estimations.size();
  
  //calculate the squared root
  rmse = rmse.array().sqrt();
  
  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  MatrixXd Hj(3,4);
  
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double px_2_plus_py_2 = fmax(pow(px, 2) + pow(py, 2), 0.0000001);
  double sqrt_px_2_plus_py_2 = sqrt(px_2_plus_py_2);
  double px_2_plus_py_2_3_2 = sqrt_px_2_plus_py_2 * px_2_plus_py_2;

  Hj << px / sqrt_px_2_plus_py_2, py / sqrt_px_2_plus_py_2, 0, 0,
        -py / px_2_plus_py_2, px / px_2_plus_py_2, 0, 0,
        py * (vx * py - vy * px) / px_2_plus_py_2_3_2, px * (vy * px - vx * py) / px_2_plus_py_2_3_2, px / sqrt_px_2_plus_py_2, py / sqrt_px_2_plus_py_2;
  
  return Hj;
}

MatrixXd Tools::ConvertCartesianToPolar(const VectorXd& x_state) {
  MatrixXd h_x(3, 1);
  
  //recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);
  
  double px_2_plus_py_2 = fmax(pow(px, 2) + pow(py, 2), 0.0000001);
  double sqrt_px_2_plus_py_2 = sqrt(px_2_plus_py_2);
  
  h_x << sqrt_px_2_plus_py_2,
         atan2(py, px),
         (px * vx + py * vy) / sqrt_px_2_plus_py_2;
  
  return h_x;
}

void Tools::NormalizePhi(Eigen::VectorXd& y) {
  double phi = y(1);
  
  phi = atan2(sin(phi), cos(phi));
  
  y(1) = phi;
}
