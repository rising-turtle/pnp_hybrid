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
#include <Eigen/Dense>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

using namespace Eigen; 
using namespace std; 
using namespace cv; 

void run_opt_3d_2d(); 
void run_opt_hybrid(); 
double  print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 


int main(int argc, char* argv[])
{
	run_opt_hybrid(); 
    // run_opt_3d_2d(); 

    return 0; 
}

void run_opt_hybrid()
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
    
    Vector3d tij(0.05, 0.05, 0.05); 
    // Vector3d tij(0.2, 0.05, 0.3); 

    Matrix3d R = Rij.transpose(); 
    Vector3d t = -Rij.transpose()*tij; 

    vector<pair<Vector3d, Vector3d>> corrs = sim.find_corrs(R, t);
           
    // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 
    Matrix3d Rij_e; //, Rij_e_t;
    Vector3d tij_e; // tij_e_t; 

   	Matrix3d R2; //, Rij_e_t;
	Vector3d t2; // tij_e_t; 

	int cnt_2d = 400; 
    int cnt_3d = 10; 
	OptSolver opt_solver; 

    // 3d-2d 
    cv::Mat rvec, tvec; 

    int trailNum = 100;
    cv::Mat RE_3d2d= Mat::zeros(trailNum, 1,CV_64FC1);
    cv::Mat RE_hybrid= Mat::zeros(trailNum, 1,CV_64FC1);
    for(int i=0; i<trailNum-1; i++){

    	cout<<"times: "<<i+1<<endl; 
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
    	RE_3d2d.at<double>(i) = print_err("3D_2D iterative: ", Rij_e, tij_e, Rij, tij); 
		    	 
    	// solve Rij_e using the specified number of 2d features 
    	me.solvePNP_2D_2D(in_2d, Rij_e, tij_e);
    	// solve tij_e     
    	// st.solveTCeres(in_3d, Rij_e, tij_e); 
    	// optimize the [Rij_e] and [tij_e]
    	RE_hybrid.at<double>(i) = print_err("initial hybrid: ", Rij_e, tij_e, Rij, tij);
    	// opt_solver.solveCeres(in_3d, Rij_e, tij_e); 
    	// print_err("hybrid after opt: ", Rij_e, tij_e, Rij, tij);

    	cout<<endl<<endl;
	}

    cv::Scalar mean, stddev;
    
    cout << trailNum <<" times "<<cnt_3d <<" depth correspondences from "<<cnt_2d << " visual features"<<endl;

    cv::meanStdDev(RE_3d2d, mean, stddev);
    cout << "3D_2D iterative:       mean: " << mean.val[0]    <<" std: " << stddev.val[0] <<endl;
    cv::meanStdDev(RE_hybrid, mean, stddev);
    cout << "hybrid_rot iterative:  mean: " << mean.val[0]    <<" std: " << stddev.val[0] <<endl;
    return ; 
}

void run_opt_3d_2d()
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

    Matrix3d Rij_e; //, Rij_e_t;
    Vector3d tij_e; // tij_e_t; 
    cv::Mat mask; 
   	Matrix3d R2; //, Rij_e_t;
	Vector3d t2; // tij_e_t; 

    int cnt_3d = 20; 
	OptSolver opt_solver; 

    // 3d-2d 
    cv::Mat rvec, tvec; 

    for(int i=0; i<10; i++){

    	cout<<"times: "<<i<<endl; 
    	vector<pair<Vector3d, Vector3d>> corrs_noise = sim.add_noise(corrs); 

	    me.solveRelativeRT_PNP(corrs_noise, Rij_e, tij_e, mask, &rvec, &tvec);
	    print_err("3D_2D EPNP: ", Rij_e, tij_e, Rij, tij); 

	    vector<pair<Vector3d, Vector3d>> inliers_3d = getInliersIndex(corrs_noise, mask); 
	    vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers_3d, cnt_3d); 

	    // compute 3d-2d transformation 
	    // me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 

	    me.solvePNP_3D_2D_given_rt(in_3d, R2, t2, rvec, tvec); 
	    print_err("3D_2D iterative: ", R2, t2, Rij, tij); 

    	opt_solver.solveCeres(in_3d, Rij_e, tij_e); 
    	print_err("Opt: ",  Rij_e, tij_e, Rij, tij);

    	cout<<endl<<endl;
	}
    return ; 

}


double print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	// cout<<pre<<" trans err: "<<dt.norm()<<" ros err: "<<computeAngle(dR)<<endl; 
    cout<<pre<<" ros err: "<<computeAngle(dR)<<endl; 
    return computeAngle(dR);
}
