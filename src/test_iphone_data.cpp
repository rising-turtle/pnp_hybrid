/*
	Feb. 8, 2021, He Zhang, fuyinzh@gmail.com 

	test iphone12's data for pnp 

*/

#include "sim_corres.h"
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

string fname1("../../debug_files/pts_image1.log"); 
string fname2("../../debug_files/pts_image2.log"); 

double FX = 1460.2; 
double FY = 1460.2; 
double CX = 965.24; 
double CY = 653.62; 

vector<Vector3d> readfile(string fname); 

void run_with_iphone_data(); 

int main(int argc, char* argv[])
{
	if(argc >= 3){
		fname1 = argv[1]; 
		fname2 = argv[2]; 
	}
	cout <<" fname1: "<<fname1<<endl; 
	cout <<" fname2: "<<fname2<<endl; 

	run_with_iphone_data(); 

	return 0; 
}


void run_with_iphone_data(){

	vector<Vector3d> pts_i = readfile(fname1); 
	vector<Vector3d> pts_j = readfile(fname2); 

	vector<pair<Vector3d, Vector3d>> corrs; 
	vector<double> vdx, vdy; 
    for(int i=0; i<pts_i.size(); i++){
        Vector3d pti = pts_i[i]; 
        Vector3d ptj = pts_j[i];
        // pti(2) += 30.; 
       	ptj(2) = 1.; 
        corrs.push_back(make_pair(pti, ptj));  
        if(i < 10)
        	cout<<"pti: "<<pti.transpose()<<" ptj: "<<ptj.transpose()<<endl; 
        vdx.push_back(fabs(pti(0) - ptj(0))); 
        vdy.push_back(fabs(pti(1) - ptj(1)));    
    }

    pair<double, double> sta_dx = getMeanStd(vdx); 
    pair<double, double> sta_dy = getMeanStd(vdy); 
    cout<<fixed<<" du: mean: "<<sta_dx.first<<" std: "<<sta_dx.second<<endl;
    cout << "dv: mean: "<<sta_dy.first<<" std: "<<sta_dy.second<<endl;

   // 3d-2d 
    MotionEstimator me;
    Matrix3d R12; 
    Vector3d t12; 
    // 3d-2d result 
    // me.solvePNP_3D_2D(corrs, R12, t12); 

    OptSolver opt_solver; 

	int cnt_3d = 10;
    // for(int cnt_3d = 10; cnt_3d < 200; cnt_3d += 10)
	{
	    vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs, cnt_3d); 
	    R12 << 0.87697,  0.311094, -0.366257,
	            -0.211306,  0.934174,  0.287522,
	            0.431594, -0.174756,  0.884979;
	    t12 << -0.528285, -1.08058, -0.619827;
	    opt_solver.solveCeres(in_3d, R12, t12); 
	    cout <<"cnt_3d: "<<cnt_3d<<endl; 
	    cout <<std::fixed<<" R12: "<<endl<< R12<<endl;
	    cout <<std::fixed<<" angle: "<<computeAngle(R12)<<endl; 
	    cout<<" t12: "<<t12.transpose()<<" t12 norm: "<<t12.norm()<<endl; 
	}
    return ; 
}

vector<Vector3d> readfile(string fname){

	ifstream inf(fname.c_str()); 

	vector<Vector3d> ret; 

	double u, v, d; 
	double x, y; 
	int i =0 ;
	while(inf.good()){

		stringstream ss; 
		string line; 
		getline(inf, line); 
		if(line.empty()) break; 

		// ss << line; 
		// ss >> u >> v >> d;

		sscanf(line.c_str(), "%lf,%lf,%lf\n", &u, &v, &d); 

		// cout <<++i<<" u: "<<u<<" v: "<<v<<" d: "<<d<<endl; 

		x = (u - CX)/FX; 
		y = (v - CY)/FY; 
		ret.push_back(Vector3d{x,y,d}); 
	}

	return ret; 

}