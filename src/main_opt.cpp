/*
    Jan. 9th, 2020, He Zhang, hzhang8@vcu.edu 

    simulate for the optimization method with projected factors 

*/

#include "sim_corres.h"
#include "solve_5pts.h"
#include "solve_translate.h"
#include "utility.h"
#include "opt_solver.h"
#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

using namespace Eigen; 
using namespace std; 

void run_opt(); 
void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 


int main(int argc, char* argv[])
{
	
    run_opt(); 

    return 0; 
}


void run_opt()
{
	SimCorr sim; 

    stringstream ss; 

    double roll(7.); // 7 
    double yaw(10.); // 10
    double pitch(5.);  // 5
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
           
    // 3d-2d 
    MotionEstimator me;

    Matrix3d Rij_e; //, Rij_e_t;
    Vector3d tij_e; // tij_e_t; 
    cv::Mat mask; 

    int cnt_3d = 20; 

    // 3d-2d 
    me.solveRelativeRT_PNP(corrs_noise, Rij_e, tij_e, mask);
    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 

    vector<pair<Vector3d, Vector3d>> inliers_3d = getInliersIndex(corrs_noise, mask); 
    vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers_3d, cnt_3d); 

    // compute 3d-2d transformation 
    me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
    print_err("3D_2D: ", Rij_e, tij_e, Rij, tij); 

    OptSolver opt_solver; 
    opt_solver.solveCeres(in_3d, Rij_e, tij_e); 
    print_err("Opt: ",  Rij_e, tij_e, Rij, tij);

    return ; 

}


void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" ros err: "<<computeAngle(dR)<<endl; 
}