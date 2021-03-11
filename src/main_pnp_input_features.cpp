/*
	Feb. 24, 2021, fuyinzh@gmail.com, He Zhang 

	run pnp with input from tracked features 

*/

#include "solve_5pts.h"
#include "solve_translate.h"
#include "utility.h"
#include "opt_solver.h"
#include <stdio.h>
#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>

#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>

using namespace Eigen; 
using namespace std; 
using namespace opengv;

string folder("../../result/feats_tracker_results"); 

vector<vector<pair<Vector3d, Vector3d>>> read_file(string filename); 

bool run_once(int cnt_3d, int cnt_2d, vector<pair<Vector3d, Vector3d>>& corres, Eigen::Matrix3d& Rij,
				 Eigen::Matrix<double, 3, 1>& tij, vector<double>& v_dR, vector<double>& v_dt);


int main(int argc, char* argv[]){

	if(argc >= 2)
		folder = string(argv[1]);

	vector<string> motion{"rot_roll", "rot_pitch", "rot_yaw", "trans_y", "trans_z"};  
	int CNT = 100; 

	string motion_record("rot_pitch"); // pitch is actually yaw 
	vector<double> gt_yaw_y{0, 3.13, 6.19, 9.25, 12.32, 15.32, 18.44}; // rot_pitch

	// vector<int> v_cnt{10, 20, 30, 40, 50, 60, 70}; 
	vector<int> v_cnt{ 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};

	for(int k=1; k<gt_yaw_y.size(); k++){

		// Ground truth
		Eigen::AngleAxisd rollAngle(D2R(0.), Eigen::Vector3d::UnitZ()); // Z
        Eigen::AngleAxisd pitchAngle(D2R(0.), Eigen::Vector3d::UnitX()); // X
		Eigen::AngleAxisd yawAngle(D2R(-gt_yaw_y[k]), Eigen::Vector3d::UnitY()); // rotate around Y 
		Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle; 
		Eigen::Matrix3d Rij = q.matrix();
        Vector3d tij = Vector3d::Zero(); 

        // 
        vector<vector<double> > vv_dR(v_cnt.size());
        for(int n=0; n<CNT; n++){

			string filename = folder + "/" + motion_record + "_" + to_string(n+1); 
			// cout<<"read file from: "<<filename<<endl; 
			vector<vector<pair<Vector3d, Vector3d>>> matched_feats = read_file(filename); 

			for(int i=0; i<v_cnt.size(); i++){
				vector<double> v_dR, v_dt; 
				// cout <<"run with features: "<<matched_feats[k].size()<<endl;
				run_once(10, v_cnt[i], matched_feats[k-1], Rij, tij, v_dR, v_dt); 
				if(v_dR[0] < 10)
					vv_dR[i].push_back(v_dR[0]); // only compute rotation error
			}
		}

		// save results 
		string out_file = motion_record + "_" + to_string(k*3) + "_hybrid.txt"; 
		ofstream ouf(out_file.c_str()); 
		for(int i=0; i<v_cnt.size(); i++){
			pair<double, double> trans = getMeanStd(vv_dR[i]);
			ouf << v_cnt[i] <<" "<< trans.first<<" "<<trans.second<<endl; 
			cout << v_cnt[i] <<" dR mean: "<< trans.first<<" std: "<<trans.second<<endl; 
		}
		ouf.close(); 
	} 
	
	return 0; 
}

bool run_once(int cnt_3d, int cnt_2d, vector<pair<Vector3d, Vector3d>>& corres, Eigen::Matrix3d& Rij,
				 Eigen::Matrix<double, 3, 1>& tij, vector<double>& v_dR, vector<double>& v_dt)
{
    // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 

    Matrix3d Rij_e, Rji_e; //, Rij_e_t;
    Vector3d tij_e; 

    vector<pair<Vector3d, Vector3d>> in_2d = getN(corres, cnt_2d); 
    vector<pair<Vector3d, Vector3d>> in_3d( in_2d.begin(), in_2d.begin()+cnt_3d);

    // ofstream ouf("tmp.log"); 
    // for(int i=0; i<in_2d.size(); i++){
    // 	ouf<<in_2d[i].first.transpose()<<" "<<in_2d[i].second.transpose()<<endl; 
    // }
    // ouf.close();

    // 3d-2d result 
    me.solvePNP_3D_2D(in_3d, Rij_e, tij_e);     
    Matrix3d dR = Rij.transpose()*Rij_e;

    // cout<<"EPNP Rij_e: "<<endl<<Rij_e<<endl; 
    // cout<<"EPNP dr: "<<R2D(computeAngle(dR))<<endl; 

    // run Hybrid 
    // 3d-2d result 
    // me.solvePNP_2D_2D(in_2d, Rij_e, tij_e); 

    // cout<<"after solvePNP_2D_2D "<<endl; 
     //derive correspondences based on random point-cloud
    bearingVectors_t points;
    bearingVectors_t bearingVectors1; 
    bearingVectors_t bearingVectors2;

    for(int j=0; j<in_2d.size(); j++){
        Vector3d pi = in_2d[j].first;
        pi(0) = pi(0)*pi(2); 
        pi(1) = pi(1)*pi(2);  
        Vector3d pj = in_2d[j].second; 
        // bearingVectors1.push_back(pi/pi.norm()); 
        points.push_back(pi); 
        bearingVectors1.push_back(pi/pi.norm());
        bearingVectors2.push_back(pj/pj.norm()); 
    }

    relative_pose::CentralRelativeAdapter adapter_rbs(
          bearingVectors1,
          bearingVectors2,
          Rij_e ); // rotation);

    // first use eight pts to compute initial rotation 
    Rij_e =  relative_pose::eigensolver(adapter_rbs);       
    dR = Rij.transpose()*Rij_e;    
    // cout<<"Rij: "<<endl<<Rij<<endl; 
    // cout<<"Rij_e: "<<endl<<Rij_e<<endl; 
    cout<<"dr: "<<R2D(computeAngle(dR))<<endl; 
    v_dR.push_back(R2D(computeAngle(dR)));  
    return true; 
}


vector<vector<pair<Vector3d, Vector3d>>> read_file(string filename)
{

	ifstream inf(filename.c_str()); 

	double depth; 
	double px, py; 

	vector<vector<pair<Vector3d, Vector3d>>> ret(6); // 6 matches, 
	while(inf.good()){

		stringstream ss; 
		string line; 
		getline(inf, line); 
		if(line.empty()) break; 
		ss << line; 
		ss >> depth >> px >> py; 

		Vector3d pti{px, py, depth}; 
		for(int j=0; j<6; j++){
			ss >> px >> py; 
			Vector3d ptj{px, py, 1.}; 
			ret[j].push_back(make_pair(pti, ptj)); 
		}
	}
	return ret; 
}

