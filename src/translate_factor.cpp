/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu 

	factors to solve t in [R,t] given R using ceres 
	
*/

#include "translate_factor.h"

TranslateFactor::TranslateFactor(const Eigen::Matrix3d& R, const Eigen::Vector3d& xi, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi):
pt_xi(xi), pt_xnj(xnj), Rji(R), cov_pt_xi(cov_xi)
{}

// bool TranslateFactor::error(Vector3d& tv, Vector2d& res)
// {

// 	Eigen::Matrix<double, 1, 3> R1 = Rji.row(0); 
// 	Eigen::Matrix<double, 1, 3> R2 = Rji.row(1); 
// 	Eigen::Matrix<double, 1, 3> R3 = Rji.row(2); 

// 	Eigen::Matrix<double, 2, 3> A; 
// 	A << 1, 0, -pt_xnj(0),
// 		 0, 1, -pt_xnj(1); 
// 	Eigen::Matrix<double, 2, 1> b; 
// 	Eigen::Matrix<double, 1, 3> dr_dx = (R1 - pt_xnj(0)*R3);
// 	Eigen::Matrix<double, 1, 3> dr_dy = (R2 - pt_xnj(1)*R3); 
// 	b(0) = (pt_xnj(0)*R3 - R1)*pt_xi; 
// 	b(1) = (pt_xnj(1)*R3 - R2)*pt_xi;

// 	// weight 
// 	Eigen::Matrix2d W = Eigen::Matrix2d::Identity(); 
// 	W(0,0) = sqrt(dr_dx * cov_pt_xi * dr_dx.transpose()); 
// 	W(1,1) = sqrt(dr_dy * cov_pt_xi * dr_dy.transpose()); 

// 	res = A*tv - b; 
// 	res = W * res;
// 	return true; 
// }

bool TranslateFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	Eigen::Vector3d T(parameters[0][0], parameters[0][1], parameters[0][2]); 
	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals); 

	Eigen::Matrix<double, 1, 3> R1 = Rji.row(0); 
	Eigen::Matrix<double, 1, 3> R2 = Rji.row(1); 
	Eigen::Matrix<double, 1, 3> R3 = Rji.row(2); 

	Eigen::Matrix<double, 2, 3> A; 
	A << 1, 0, -pt_xnj(0),
		 0, 1, -pt_xnj(1); 
	Eigen::Matrix<double, 2, 1> b; 
	Eigen::Matrix<double, 1, 3> dr_dx = (R1 - pt_xnj(0)*R3);
	Eigen::Matrix<double, 1, 3> dr_dy = (R2 - pt_xnj(1)*R3); 
	b(0) = (pt_xnj(0)*R3 - R1)*pt_xi; 
	b(1) = (pt_xnj(1)*R3 - R2)*pt_xi;

	// weight 
	Eigen::Matrix2d W = Eigen::Matrix2d::Identity(); 
	W(0,0) = sqrt(dr_dx * cov_pt_xi * dr_dx.transpose()); 
	W(1,1) = sqrt(dr_dy * cov_pt_xi * dr_dy.transpose()); 

	residual = A*T - b; 
	residual = W * residual;

	if(jacobians){
		if(jacobians[0]){
			Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_t(jacobians[0]); 

			jacobian_t = W*A;

		}
	}

	return true; 
}