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
    	residuals[0] =100.* (xp - T(observed_u));
    	residuals[1] =100.* (yp - T(observed_v));
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

struct SampsonEssential
{
	SampsonEssential(double ui, double vi, double uj, double vj)
		:m_ui(ui), m_vi(vi), m_uj(uj), m_vj(vj)
		{}

	template <typename T>
	bool operator()(const T* const camera_quaternion, const T* const camera_Tji, T* residuals) const
	{
		T R[9]; // R is Rji camera_quaternion is Qji 
		ceres::QuaternionToRotation(camera_quaternion, R); 

		T camera_T[3]; // tij = -Rji'*tji

		camera_T[0] = -R[0]*camera_Tji[0] - R[3]*camera_Tji[1] - R[6]*camera_Tji[2]; 
		camera_T[1] = -R[1]*camera_Tji[0] - R[4]*camera_Tji[1] - R[7]*camera_Tji[2]; 
		camera_T[2] = -R[2]*camera_Tji[0] - R[5]*camera_Tji[1] - R[8]*camera_Tji[2]; 


		T t[9]; // t = [tij]
		t[0] = T(0); t[1] = -camera_T[2]; t[2] = camera_T[1]; 
		t[3] = camera_T[2]; t[4] = T(0); t[5] = -camera_T[0]; 
		t[6] = -camera_T[1]; t[7] = camera_T[0]; t[8] = T(0); 
		T E[9]; // E = Rji * [tij]
		E[0] = R[0]*t[0] + R[1]*t[3] + R[2]*t[6]; 
		E[1] = R[0]*t[1] + R[1]*t[4] + R[2]*t[7];
		E[2] = R[0]*t[2] + R[1]*t[5] + R[2]*t[8];

		E[3] = R[3]*t[0] + R[4]*t[3] + R[5]*t[6]; 
		E[4] = R[3]*t[1] + R[4]*t[4] + R[5]*t[7];
		E[5] = R[3]*t[2] + R[4]*t[5] + R[5]*t[8];

		E[6] = R[6]*t[0] + R[7]*t[3] + R[8]*t[6]; 
		E[7] = R[6]*t[1] + R[7]*t[4] + R[8]*t[7];
		E[8] = R[6]*t[2] + R[7]*t[5] + R[8]*t[8];

		T ep_i[3]; // ep_i = E * [ui, vi, 1]
		ep_i[0] = E[0]*T(m_ui) + E[1]*T(m_vi) + E[2]; 
		ep_i[1] = E[3]*T(m_ui) + E[4]*T(m_vi) + E[5]; 
		ep_i[2] = E[6]*T(m_ui) + E[7]*T(m_vi) + E[8];

		T ep_j[3]; // ep_j = E' * [uj, vj, 1]
		ep_j[0] = E[0]*T(m_uj) + E[3]*T(m_vj) + E[6]; 
		ep_j[1] = E[1]*T(m_uj) + E[4]*T(m_vj) + E[7]; 
		ep_j[2] = E[2]*T(m_uj) + E[5]*T(m_vj) + E[8];

		T epsilon; // epsilon = [uj, vj, 1]' E [ui, vi, 1]
		epsilon = T(m_uj) * ep_i[0] + T(m_vj) * ep_i[1] + ep_i[2]; 

		T JJ; // J'J 
		JJ = ep_i[0]*ep_i[0] + ep_i[1]*ep_i[1] + ep_j[0]*ep_j[0] + ep_j[1]*ep_j[1]; 
		T inv_JJ = 1./JJ; 
		residuals[0] = - ep_j[0]*inv_JJ*epsilon; 
		residuals[1] = - ep_j[1]*inv_JJ*epsilon;
		residuals[2] = - ep_i[0]*inv_JJ*epsilon;
		residuals[3] = - ep_i[1]*inv_JJ*epsilon; 

    	return true;
	}

	static ceres::CostFunction* Create(const double ui, const double vi,
	                                   const double uj, const double vj) 
	{
	  return (new ceres::AutoDiffCostFunction<
	          SampsonEssential, 4, 4, 3>(
	          	new SampsonEssential(ui, vi, uj, vj)));
	}


	double m_ui, m_vi;
	double m_uj, m_vj;
};

class OptSolver 
{
public:
	OptSolver();

	// notice here corres must be inliers 
	bool solveCeres(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &Rij, Vector3d &tv);

	// combine 2D and 3D 
	bool solveCeresHybrid(const vector<pair<Vector3d, Vector3d>> &corres_3d, 
		const vector<pair<Vector3d, Vector3d>> &corres_2d, Matrix3d &Rij, Vector3d &tv);

	// camera pose 
	double m_tji[1][3]; 
	double m_qji[1][4]; 

	double para_pt[1000][3]; 

	int feature_num;
};