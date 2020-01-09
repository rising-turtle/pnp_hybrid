/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu 

	solve t in [R,t] given R 
	
*/

#pragma once


#include <ceres/ceres.h>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;  

class SolveTranslate 
{
public:
	
	SolveTranslate(); 

	// notice here corres must be inliers 
	bool solveTCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tv);
	bool solveTProjCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tij);

	bool solveTScaleCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, const Vector3d& ntij, Vector3d &tv);
	bool solveTProjScaleCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, const Vector3d& ntij, Vector3d &tij);

	bool solveTCeresWithPt(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tij);
	bool solveTProjCeresWithPt(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tij);

	double para_T[1][3]; // for translation 
	double para_s[1][0]; // for scale 
	double para_pt[1000][3]; 

};