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
#include "utility.h"

using namespace std; 
using namespace Eigen; 

class SimCorr
{
public:
    SimCorr(); 
    ~SimCorr(); 
    
    // init feature points randomly locates in range z-axis [2-7] meters 
    map<int, Vector3d > mg_feats; // globale feature points [normalized_u, normalized_v, z]
    void init(int pix_step = 10); // first pose is at [0, 0, 0] towards z axis 
    
    // add noise to the corresponse features 
    vector<pair<Vector3d, Vector3d>> add_noise( vector<pair<Vector3d, Vector3d>>& corres, double pix_std = 0.5/460., poly d_std=poly());

    // find correspond of the 
    vector<pair<Vector3d, Vector3d>> find_corrs(Matrix3d& R, Vector3d& t); 

    // follow opengv to generate random points 3D-2D style 
    vector<pair<Vector3d, Vector3d>> addNoise3D2D( vector<pair<Vector3d, Vector3d>>& corres, double pix_std = 0.5/460., poly d_std=poly()); 

    // follow opengv to add noise in the normal plane of the bearingVector 
    Eigen::Vector3d addNoiseRay( double noiseLevel, Eigen::Vector3d cleanPoint ); 
};

