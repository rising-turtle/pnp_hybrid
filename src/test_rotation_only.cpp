/*
	Jan. 28, 2021, He Zhang, fuyinzh@gmail.com
	
	test relative rotation given two sets of points  

*/

#include "rotation_only.h"

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

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 
void run_test();

int main(int argc, char* argv[])
{
	run_test(); 
    return 0; 
}

void run_test()
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

    // cout<<"ground truth Rij: "<<endl<<Rij<<endl; 
    
    Vector3d tij(0.2, 0.05, 0.3); 

    Matrix3d R = Rij.transpose(); 
    Vector3d t = -Rij.transpose()*tij; 

    vector<pair<Vector3d, Vector3d>> corrs = sim.find_corrs(R, t);
     // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 
    Matrix3d Rij_e; //, Rij_e_t;
    Vector3d tij_e; // tij_e_t; 

   	// rotation only 
   	RotationOnly ro; 
   	Matrix3d Rij_o; 

	int cnt_2d = 40; 
    int cnt_3d = 8; //20; 

    // 3d-2d 
    cv::Mat rvec, tvec; 

    for(int i=0; i<10; i++){

    	cout<<"times: "<<i<<endl; 
    	vector<pair<Vector3d, Vector3d>> corrs_noise = sim.add_noise(corrs); 

    	// 2d-2d get rotation 
    	cv::Mat mask, mask_t; 
    	me.solveRelativeRT(corrs_noise, Rij_e, tij_e, mask); 

    	// with the speficied number of features 
    	vector<pair<Vector3d, Vector3d>> inliers = getInliers(corrs_noise, mask); 

    	// run 3d-2d find out these 
    	me.solveRelativeRT_PNP(inliers, Rij_e, tij_e, mask_t);
    	vector<pair<Vector3d, Vector3d>> inliers_3d = getInliersIndex(inliers, mask_t); 
    	vector<pair<Vector3d, Vector3d>> in_2d = getN(inliers, cnt_2d); 
    	vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers_3d, cnt_3d); 

    	// 3d-2d result 
    	me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
    	print_err("EPNP : ", Rij_e, tij_e, Rij, tij); 
		cout<<"EPNP: Rij: "<<endl<<Rij_e<<endl; 

		ro.computeRotationSac(in_2d, Rij_o); 
		cout<<"Rotation only: Rij: "<<endl<<Rij_o<<endl; 
		print_err("Rotation only: ", Rij_o, tij_e, Rij, tij); 

    	cout<<endl<<endl;
	}
    return ; 

}


void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" ros err: "<<computeAngle(dR)<<endl; 
}