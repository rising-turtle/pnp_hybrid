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

#define SQ(x) ((x)*(x))

#define D2R(d) (((d)/180.)*M_PI)
#define R2D(r) (((r)/M_PI)*180.)

extern vector<pair<Vector3d, Vector3d>> getN(const vector<pair<Vector3d, Vector3d>>& in, int N);
extern double computeAngle(Matrix3d& R); 

extern pair<double, double> getMeanStd(vector<double>& data); 

extern vector<pair<Vector3d, Vector3d>> getInliers(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask);     
extern vector<pair<Vector3d, Vector3d>> getInliersIndex(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask);     

extern double sum_error(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d& Rij, const Vector3d& tij); 
extern double sum_error_2d(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d& Rij, const Vector3d& tij);

struct poly{
	poly(double para[3]){
		a1 = para[0]; a2 = para[1]; a3 = para[2]; 
	}
	poly(){
	    // default struct_core depth's variance w.r.t. depth  
	    a1 = 0.00155816; a2 = -0.00362021; a3 = 0.00452812;
	    // default depth's variance fx = 525  
	    // a1 = 1./(800.*0.1); a2 = 0; a3 = 0; 
	    a2 = a3 = 0;
	}
	double y(double x){
		if(x <= 0.75)
			return 0.0007;
		return (a1*x*x + a2*x+a3);
	}
	double a1,a2,a3;
	int r,c; 
};
