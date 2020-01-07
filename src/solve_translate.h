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

	double para_T[0][3]; 

};