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
  for(int i=0; i < estimations.size(); ++i){
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
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float px_2_plus_py_2 = pow(px,2) + pow(py,2);
  float sqrt_px_2_plus_py_2 = sqrt(pow(px,2) + pow(py,2));
  //check division by zero
  if (px_2_plus_py_2 > 0.0001)
  {
    //compute the Jacobian matrix
    Hj << px / sqrt_px_2_plus_py_2, py / sqrt_px_2_plus_py_2, 0, 0,
          -py/px_2_plus_py_2, px/px_2_plus_py_2, 0, 0,
          (py*(vx*py-vy*px))/sqrt(pow(px_2_plus_py_2, 3)), (px*(vx*py-vy*px))/pow(px_2_plus_py_2, 3/2), px / sqrt_px_2_plus_py_2, py / sqrt_px_2_plus_py_2
    ;
  }
  else
  {
    std::cout << "CalculateJacobian () - Error - Dividion by Zero" << std::endl;
  }
  
  return Hj;
}

MatrixXd Tools::ConvertCartesianToPolar(const VectorXd& x_state) {
  MatrixXd h_x(3, 1);
  
  //recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float px_2_plus_py_2 = pow(px,2) + pow(py,2);
  float sqrt_px_2_plus_py_2 = sqrt(pow(px,2) + pow(py,2));
  
  //check division by zero
  if (px_2_plus_py_2 > 0.0001)
  {
    //compute the Jacobian matrix
    h_x << sqrt_px_2_plus_py_2,
           atan(py/px),
           (px * vx + py * vy) / sqrt_px_2_plus_py_2;
  }
  else
  {
    std::cout << "CalculateCartesianToPolarMatrix () - Error - Dividion by Zero" << std::endl;
  }
  
  return h_x;
}

void Tools::NormalizePhi(Eigen::VectorXd& y) {
  float phi = y(1);
  
  phi = phi - static_cast<int>(phi / (2 * M_PI));
  phi = phi > M_PI ? phi - 2 * M_PI : phi;
  
  y(1) = phi;
}
