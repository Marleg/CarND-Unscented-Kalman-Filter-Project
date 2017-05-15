#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  if(estimations.size() == 0 || estimations.size() != ground_truth.size()){
    // std::cout << "Invalid estimation or ground truth data" << std::endl;
    return rmse;
  }

  for(int unsigned i = 0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculating the mean
  rmse = rmse / estimations.size();

  //calculating the squared root
  rmse = rmse.array().sqrt();

  // std::cout << "tools:CalculateRMSE executed" << std::endl;
  return rmse;

}
