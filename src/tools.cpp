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
  TO_DID:
    * Calculate the RMSE here.
  */
  VectorXd rmse(5);
  rmse.fill(0.0);

	// check the validity of the following inputs:
	// * the estimation vector size should not be zero
	// * the estimation vector size should equal ground truth vector size
	if (estimations.size() == 0) {
		//cout << "Estimation vector size be zero!" << endl;
		return rmse;
	};

	if (estimations.size() != ground_truth.size()) {
		//cout << "Estimation and ground truth vector sizes not sames!" << endl;
		return rmse;
	};

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){

        VectorXd residual = estimations[i] - ground_truth[i];
				residual = residual.array().square();
        rmse += residual;
	}

	//calculate the mean
	rmse /= estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}
