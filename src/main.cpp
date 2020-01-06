/*
    Jan. 6th, 2020, He Zhang, hzhang8@vcu.edu 

    entry file for simulation 

*/

#include "sim_corres.h"
#include <Eigen/Core>

using namespace Eigen; 
using namespace std; 

int main(int argc, char* argv[])
{
    SimCorr sim; 

    Matrix3d R = Matrix3d::Identity(); 
    
    Vector3d t(0.2, 0.05, 0.3); 

    vector<pair<Vector3d, Vector3d>> corrs = sim.find_corrs(R, t); 
    
    /*for(int i=0; i<corrs.size(); i++){
	
	cout<<"i: "<<i<<" pi: "<<corrs[i].first.transpose()<<" pj: "<<corrs[i].second.transpose()<<endl;  
        
    }*/

    return 0; 
}
