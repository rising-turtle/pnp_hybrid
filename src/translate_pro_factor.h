/*

	Jan. 9, 2020, He Zhang, hzhang8@vcu.edu 

	factors to solve t in [R,t] given R using ceres 
	
	minize projected error 
*/

#pragma once


#include <ceres/ceres.h>
#include <Eigen/Dense>


class TranslateProFactor : public ceres::SizedCostFunction<2, 3>
{
public:
	TranslateProFactor(const Eigen::Matrix3d& R, const Eigen::Vector3d& xi, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi);
	virtual bool Evaluate(double const* const *params, double *residuals, double**jacobians) const; 

	Eigen::Vector3d pt_xi; 
	Eigen::Vector2d pt_xnj; 
	Eigen::Matrix3d Rji; 
	Eigen::Matrix3d cov_pt_xi; 
	// Eigen::Matrix2d cov_pt_xnj; 
	
};


class TranslateProScaleFactor : public ceres::SizedCostFunction<2, 1>
{
public:
	TranslateProScaleFactor(const Eigen::Matrix3d& R, const Eigen::Vector3d& nt, const Eigen::Vector3d& xi, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi);
	virtual bool Evaluate(double const* const *params, double *residuals, double**jacobians) const; 

	Eigen::Vector3d pt_xi; 
	Eigen::Vector2d pt_xnj; 
	Eigen::Matrix3d Rji; 
	Eigen::Vector3d ntji; // normalized tji 
	Eigen::Matrix3d cov_pt_xi; 
};

class TranslateProWithPtFactor : public ceres::SizedCostFunction<2, 3, 3>
{
public:
	TranslateProWithPtFactor(const Eigen::Matrix3d& R, const Eigen::Vector2d& xnj, 
			const Eigen::Matrix3d& cov_xi);
	virtual bool Evaluate(double const* const *params, double *residuals, double**jacobians) const; 

	Eigen::Vector2d pt_xnj; 
	Eigen::Matrix3d Rji; 
	Eigen::Matrix3d cov_pt_xi; 
};