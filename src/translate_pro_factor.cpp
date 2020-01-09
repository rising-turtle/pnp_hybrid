/*

	Jan. 9, 2020, He Zhang, hzhang8@vcu.edu 

	factors to solve t in [R,t] given R using ceres 

	with projected factors 
	
*/

#include "translate_pro_factor.h"
#include "utility.h"
using namespace Eigen; 

TranslateProWithPtFactor::TranslateProWithPtFactor(const Eigen::Matrix3d& R, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi):
			pt_xnj(xnj), Rji(R), cov_pt_xi(cov_xi)
{}

bool TranslateProWithPtFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	Eigen::Vector3d T(parameters[0][0], parameters[0][1], parameters[0][2]); 
	Eigen::Vector3d pt_xi(parameters[1][0], parameters[1][1], parameters[1][2]); 
	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals); 

	Vector3d pt_j = Rji*pt_xi + T; 
	Vector2d pt_nj(pt_j.x()/pt_j.z(), pt_j.y()/pt_j.z()); 
	residual = pt_nj - pt_xnj; 

	Eigen::Matrix<double, 1, 3> R1 = Rji.row(0); 
	Eigen::Matrix<double, 1, 3> R2 = Rji.row(1); 
	Eigen::Matrix<double, 1, 3> R3 = Rji.row(2); 

	// weight 
	Eigen::Matrix2d W = Eigen::Matrix2d::Identity();

	Eigen::Matrix<double, 1, 3> dex_dpi, dey_dpi; 
	dex_dpi = R1/pt_j.z() - (pt_j.x()/(SQ(pt_j.z())))*R3;
	dey_dpi = R2/pt_j.z() - (pt_j.y()/(SQ(pt_j.z())))*R3;

	// W(0,0) = 1./sqrt(dex_dpi * cov_pt_xi * dex_dpi.transpose()); 
	// W(1,1) = 1./sqrt(dey_dpi * cov_pt_xi * dey_dpi.transpose()); 

	double sum_trace = cov_pt_xi(0,0) + cov_pt_xi(1,1) + cov_pt_xi(2,2); 
	W(0,0) = 1./sum_trace; 
	W(1,1) = 1./sum_trace; 

	residual = W * residual;

	if(jacobians){
		if(jacobians[0]){
			Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_t(jacobians[0]); 

			jacobian_t << 1./pt_j.z(), 0, -pt_j.x()/SQ(pt_j.z()), 
							0, 1./pt_j.z(), -pt_j.y()/SQ(pt_j.z()); 
			jacobian_t = W * jacobian_t;
		}
		if(jacobians[1]){
			Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_t(jacobians[1]);

			jacobian_t.row(0) = dex_dpi; 
			jacobian_t.row(1) = dey_dpi; 
			jacobian_t = W * jacobian_t; 
		}
	}

	return true; 
}

TranslateProFactor::TranslateProFactor(const Eigen::Matrix3d& R, const Eigen::Vector3d& xi, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi):
pt_xi(xi), pt_xnj(xnj), Rji(R), cov_pt_xi(cov_xi)
{}

bool TranslateProFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	Eigen::Vector3d T(parameters[0][0], parameters[0][1], parameters[0][2]); 
	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals); 

	Vector3d pt_j = Rji*pt_xi + T; 
	Vector2d pt_nj(pt_j.x()/pt_j.z(), pt_j.y()/pt_j.z()); 
	residual = pt_nj - pt_xnj; 

	// weight 
	Eigen::Matrix2d W = Eigen::Matrix2d::Identity(); 

	double sum_trace = cov_pt_xi(0,0) + cov_pt_xi(1,1) + cov_pt_xi(2,2); 
	W(0,0) = 1./sum_trace; 
	W(1,1) = 1./sum_trace; 
	// Eigen::Matrix<double, 1, 3> R1 = Rji.row(0); 
	// Eigen::Matrix<double, 1, 3> R2 = Rji.row(1); 
	// Eigen::Matrix<double, 1, 3> R3 = Rji.row(2); 

	// Eigen::Matrix<double, 1, 3> dex_dpi, dey_dpi; 
	// dex_dpi = R1/pt_j.z() - (pt_j.x()/(SQ(pt_j.z())))*R3;
	// dey_dpi = R2/pt_j.z() - (pt_j.y()/(SQ(pt_j.z())))*R3;

	// W(0,0) = 1./sqrt(dex_dpi * cov_pt_xi * dex_dpi.transpose()); 
	// W(1,1) = 1./sqrt(dey_dpi * cov_pt_xi * dey_dpi.transpose()); 

	residual = W * residual;

	if(jacobians){
		if(jacobians[0]){
			Eigen::Map<Eigen::Matrix<double, 2, 3, Eigen::RowMajor>> jacobian_t(jacobians[0]); 

			jacobian_t << 1./pt_j.z(), 0, -pt_j.x()/SQ(pt_j.z()), 
							0, 1./pt_j.z(), -pt_j.y()/SQ(pt_j.z()); 
			jacobian_t = W * jacobian_t;
		}
	}

	return true; 
}


TranslateProScaleFactor::TranslateProScaleFactor(const Eigen::Matrix3d& R, const Eigen::Vector3d& nt, const Eigen::Vector3d& xi, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi):
pt_xi(xi), pt_xnj(xnj), Rji(R), ntji(nt), cov_pt_xi(cov_xi)
{}

bool TranslateProScaleFactor::Evaluate(double const *const *parameters, double *residuals, double **jacobians) const
{
	double scale = parameters[0][0]; 
	Eigen::Vector3d T = scale * ntji; 
	Eigen::Map<Eigen::Matrix<double, 2, 1>> residual(residuals); 

	Vector3d pt_j = Rji*pt_xi + T; 
	Vector2d pt_nj(pt_j.x()/pt_j.z(), pt_j.y()/pt_j.z()); 
	residual = pt_nj - pt_xnj; 

	// weight 
	Eigen::Matrix2d W = Eigen::Matrix2d::Identity(); 

	// double sum_trace = cov_pt_xi(0,0) + cov_pt_xi(1,1) + cov_pt_xi(2,2); 
	// W(0,0) = 1./sum_trace; 
	// W(1,1) = 1./sum_trace; 
	Eigen::Matrix<double, 1, 3> R1 = Rji.row(0); 
	Eigen::Matrix<double, 1, 3> R2 = Rji.row(1); 
	Eigen::Matrix<double, 1, 3> R3 = Rji.row(2); 
	Eigen::Matrix<double, 1, 3> dex_dpi, dey_dpi; 
	dex_dpi = R1/pt_j.z() - (pt_j.x()/(SQ(pt_j.z())))*R3;
	dey_dpi = R2/pt_j.z() - (pt_j.y()/(SQ(pt_j.z())))*R3;

	W(0,0) = 1./sqrt(dex_dpi * cov_pt_xi * dex_dpi.transpose()); 
	W(1,1) = 1./sqrt(dey_dpi * cov_pt_xi * dey_dpi.transpose()); 

	residual = W * residual;

	if(jacobians){
		if(jacobians[0]){
			Eigen::Map<Eigen::Matrix<double, 2, 1>> jacobian_t(jacobians[0]); 

			Eigen::Matrix<double, 2, 3> de_dT;
			de_dT << 1./pt_j.z(), 0, -pt_j.x()/SQ(pt_j.z()), 
							0, 1./pt_j.z(), -pt_j.y()/SQ(pt_j.z()); 
			jacobian_t = W * de_dT * ntji;
		}
	}

	return true; 
}