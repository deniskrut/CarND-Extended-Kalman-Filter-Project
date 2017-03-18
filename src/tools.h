#ifndef TOOLS_H_
#define TOOLS_H_
#include <vector>
#include "Eigen/Dense"

class Tools {
public:
  /**
  * Constructor.
  */
  Tools();

  /**
  * Destructor.
  */
  virtual ~Tools();

  /**
  * A helper method to calculate RMSE.
  */
  Eigen::VectorXd CalculateRMSE(const std::vector<Eigen::VectorXd> &estimations, const std::vector<Eigen::VectorXd> &ground_truth);

  /**
  * A helper method to calculate Jacobians.
  */
  Eigen::MatrixXd CalculateJacobian(const Eigen::VectorXd& x_state);
  
  /**
   * A helper method to calculate Cartesian to Polar conversion matrix.
   */
  Eigen::MatrixXd CalculateCartesianToPolarMatrix(const Eigen::VectorXd& x_state);
  
  /**
   * A helper method to normalize ϕ.
   */
  void NormalizePhi(Eigen::VectorXd& y);

};

#endif /* TOOLS_H_ */
