/*
    Jan. 6th, 2020, He Zhang, hzhang8@vcu.edu 

    entry file for simulation 

*/

#include "sim_corres.h"
#include "solve_5pts.h"
#include <Eigen/Core>
#include <cmath>

using namespace Eigen; 
using namespace std; 

#define D2R(d) (((d)/180.)*M_PI)
#define R2D(r) (((r)/M_PI)*180.)

int main(int argc, char* argv[])
{
    SimCorr sim; 

    double roll(7.); 
    double yaw(10.);
    double pitch(5.); 
    cout<<" ground truch yaw: "<<yaw<<" pitch: "<<pitch<<" roll: "<<roll<<endl;

    Eigen::AngleAxisd rollAngle(D2R(roll), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd yawAngle(D2R(yaw), Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd pitchAngle(D2R(pitch), Eigen::Vector3d::UnitY());
    // Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle; 
    Eigen::Matrix3d Rij = q.matrix();

    cout<<"ground truth Rij: "<<endl<<Rij<<endl; 
    
    Vector3d tij(0.2, 0.05, 0.3); 

    Matrix3d R = Rij.transpose(); 
    Vector3d t = -Rij.transpose()*tij; 

    vector<pair<Vector3d, Vector3d>> corrs = sim.find_corrs(R, t); 
    vector<pair<Vector3d, Vector3d>> corrs_noise = sim.add_noise(corrs); 

    /*for(int i=0; i<corrs.size(); i++){
	
	cout<<"i: "<<i<<" pi: "<<corrs[i].first.transpose()<<" pj: "<<corrs[i].second.transpose()<<endl;  
        
    }*/
    MotionEstimator me;
    Matrix3d Rij_e;
    Vector3d tij_e;  
    me.solveRelativeRT(corrs_noise, Rij_e, tij_e); 

    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 
    cout << "to Euler angles:" << endl;
    cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    cout <<"main.cpp: 2d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 

    me.solveRelativeRT_PNP(corrs_noise, Rij_e, tij_e);
    ea = Rij_e.eulerAngles(2, 1, 0); 

    cout << "to Euler angles:" << endl;
    cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    cout <<"main.cpp: 3d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 

    return 0; 
}
