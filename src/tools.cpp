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
  if (px != 0 && py != 0)
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
