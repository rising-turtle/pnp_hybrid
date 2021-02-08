/*
    Feb. 4th, 2021, He Zhang, fuyinzh@gmail.com 

    compare opengv pnps 

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

#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>

using namespace Eigen; 
using namespace std; 
using namespace opengv; 

void run_opt_3d_2d(); 
void run_opt_hybrid(); 
void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 


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
    
    Vector3d tij(0.2, 0.05, 0.3); 

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

	int cnt_2d = 40; 
    int cnt_3d = 20; 
	OptSolver opt_solver; 

    // 3d-2d 
    cv::Mat rvec, tvec; 

    for(int i=0; i<10; i++){

    	cout<<"times: "<<i<<endl; 
        vector<pair<Vector3d, Vector3d>> corrs_noise = sim.addNoise3D2D(corrs); 
    	vector<pair<Vector3d, Vector3d>> in_2d = getN(corrs_noise, cnt_2d); 
    	vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs_noise, cnt_3d); 

    	// 3d-2d result 
    	me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
    	print_err("3D_2D iterative: ", Rij_e, tij_e, Rij, tij); 
		

        // run 2D-2D result, eigen solver 

          //derive correspondences based on random point-cloud
        bearingVectors_t bearingVectors1;
        bearingVectors_t bearingVectors2;

        for(int j=0; j<in_2d.size(); j++){
            Vector3d& pi = in_2d[j].first; 
            bearingVectors1.push_back(pi/pi.norm()); 
            bearingVectors2.push_back(in_2d[j].second); 
        }
        
        Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity(); 
        //create a central relative adapter
        relative_pose::CentralRelativeAdapter adapter(
              bearingVectors1,
              bearingVectors2,
              rotation);

        Rij_e =  relative_pose::eigensolver(adapter);  
    	st.solveTCeres(in_3d, Rij_e, tij_e); 
    	// optimize the [Rij_e] and [tij_e]
    	print_err("initial hybrid: ", Rij_e, tij_e, Rij, tij);
    	opt_solver.solveCeres(in_3d, Rij_e, tij_e); 
    	print_err("hybrid after opt: ", Rij_e, tij_e, Rij, tij);

    	cout<<endl<<endl;
	}
    return ; 
}



void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" rotation err: "<<computeAngle(dR)<<endl; 
}