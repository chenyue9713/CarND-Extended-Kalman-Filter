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
  rmse << 0,0,0,0;
  
  if(estimations.size() == 0)
  {
    cout << "The estimation vector size should not be zero!" << endl;
    return rmse;
  }
  
  if(estimations.size() != ground_truth.size())
  {
    cout << "The estimation vector size should equal to size of ground Truth vector !" << endl;
    return rmse;
  }
  
  for(unsigned int i = 0; i < estimations.size(); ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }
  
  rmse = rmse / estimations.size();
  
  rmse = rmse.array().sqrt();
  if (rmse[0] > 0.11 ||
      rmse[1] > 0.11 ||
      rmse[2] > 0.52 ||
      rmse[3] > 0.52)
  {
  	cout << "rmse exceed the limiting value" << endl;
    cout << "rmse: " << rmse[0] << ", " << rmse[1] << ", "<< rmse[2] << ", "<< rmse[3] << endl;
  }
  
  
  return rmse;
  
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  /**
  TODO:
    * Calculate a Jacobian here.
  */
  MatrixXd Hj(3,4);
  
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  
  float c1 = px * px + py * py;
  float c2 = sqrt(c1);
  float c3 = c1 * c2;
  
  Hj << 1, 1, 0, 0,
    	1, 1, 0, 0,
   		1, 1, 1, 1;
  
  if(fabs(c1) < 0.0001)
  {
    cout << "CalculateJacobian() - Error -Division by Zero" << endl;
    return Hj;
  }
  
  Hj << (px/c2), (py/c2), 0, 0,
    	-(py/c1), (px/c1),0, 0,
   		py*(vx*py-vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;
  
  return Hj;
  
}
