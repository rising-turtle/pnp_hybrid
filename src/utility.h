/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu  

	some utility functions 

*/

#pragma once 

#include <Eigen/Core>
#include <vector>
#include <map>
#include <opencv2/opencv.hpp>

using namespace std; 
using namespace Eigen;

extern vector<pair<Vector3d, Vector3d>> getN(const vector<pair<Vector3d, Vector3d>>& in, int N);
extern double computeAngle(Matrix3d& R); 

extern pair<double, double> getMeanStd(vector<double>& data); 

extern vector<pair<Vector3d, Vector3d>> getInliers(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask);     
extern vector<pair<Vector3d, Vector3d>> getInliersIndex(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask);     
