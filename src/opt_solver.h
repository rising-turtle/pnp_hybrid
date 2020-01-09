/*
	Jan. 9, 2020, He Zhang, hzhang8@vcu.edu 

	functions related to sfm 

*/

#pragma once

#include <ceres/ceres.h>
#include <ceres/rotation.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cstdlib>
#include <deque>
#include <map>
#include <opencv2/core/eigen.hpp>
#include <opencv2/opencv.hpp>
using namespace Eigen;
using namespace std;

struct SFMFeature
{
    bool state;
    int id;
    vector<pair<int,Vector2d>> observation;
    double position[3];
    double depth;
};

struct ReprojectionError3D
{
	ReprojectionError3D(double observed_u, double observed_v)
		:observed_u(observed_u), observed_v(observed_v)
		{}

	template <typename T>
	bool operator()(const T* const camera_R, const T* const camera_T, const T* point, T* residuals) const
	{
		T p[3];
		ceres::QuaternionRotatePoint(camera_R, point, p);
		p[0] += camera_T[0]; p[1] += camera_T[1]; p[2] += camera_T[2];
		T xp = p[0] / p[2];
    	T yp = p[1] / p[2];
    	residuals[0] = xp - T(observed_u);
    	residuals[1] = yp - T(observed_v);
    	return true;
	}

	static ceres::CostFunction* Create(const double observed_x,
	                                   const double observed_y) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          ReprojectionError3D, 2, 4, 3, 3>(
	          	new ReprojectionError3D(observed_x,observed_y)));
	}

	double observed_u;
	double observed_v;
};

class OptSolver 
{
public:
	OptSolver();

	// notice here corres must be inliers 
	bool solveCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tv);

	// camera pose 
	double m_tji[1][3]; 
	double m_qji[1][4]; 

	double para_pt[1000][3]; 

	int feature_num;
};