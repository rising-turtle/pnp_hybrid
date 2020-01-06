/*

    Jan. 6th, 2020, He Zhang, hzhang8@vcu.edu 
    
    simulate point correspondences given [R, t]

*/

#pragma once 

#include <iostream>
#include <map>
#include <vector>
#include <Eigen/Core>
#include <Eigen/Geometry>

using namespace std; 
using namespace Eigen; 

struct poly{
	poly(double para[3]){
		a1 = para[0]; a2 = para[1]; a3 = para[2]; 
	}
	poly(){
	    // default struct_core depth's variance w.r.t. depth  
	    a1 = 0.00155816; a2 = -0.00362021; a3 = 0.00452812;
	}
	double y(double x){
		if(x <= 0.75)
			return 0.0007;
		return (a1*x*x + a2*x+a3);
	}
	double a1,a2,a3;
	int r,c; 
};

class SimCorr
{
public:
    SimCorr(); 
    ~SimCorr(); 
    
    // init feature points randomly locates in range z-axis [2-7] meters 
    map<int, Vector3d > mg_feats; // globale feature points [normalized_u, normalized_v, z]
    void init(int pix_step = 5); // first pose is at [0, 0, 0] towards z axis 
    
    // add noise to the corresponse features 
    vector<pair<Vector3d, Vector3d>> add_noise( vector<pair<Vector3d, Vector3d>>& corres, double pix_std = 0.5/460., poly d_std=poly());

    // find correspond of the 
    vector<pair<Vector3d, Vector3d>> find_corrs(Matrix3d& R, Vector3d& t); 

};

