/*
    Feb. 7th, 2021, He Zhang, fuyinzh@gmail.com 

    monte carlo against different noise  

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
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include "random_generators.hpp"
#include "experiment_helpers.hpp"

using namespace Eigen; 
using namespace std; 
using namespace opengv; 

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 

bool run_once(double noise, int num, vector<double>& dR_epnp, vector<double>& dt_epnp); 

void run_monte_carlo( vector<double> v_noise, int cnt_points = 10,  int TIMES= 1000) ;

int main(int argc, char* argv[])
{
    vector<double> v_noise{0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0};
    run_monte_carlo(v_noise, 10, 2000);
    return 0; 
}


void run_monte_carlo( vector<double> v_noise, int cnt_points,  int TIMES)
{
    ofstream ouf_epnp("NOISE_EPNP.log"); 
    ofstream ouf_hybrid("NOISE_HYBRID.log"); 

    for(int i=0; i<v_noise.size(); i++){

        vector<double> hybrid_et, hybrid_er, epnp_et, epnp_er; 

        for(int k=0; k < TIMES; ){
            vector<double> v_dr, v_dt; 
            if(run_once(v_noise[i], cnt_points, v_dr, v_dt)){
                k++; 
                epnp_er.push_back(v_dr[0]); 
                epnp_et.push_back(v_dt[0]);
                hybrid_er.push_back(v_dr[1]);
                hybrid_et.push_back(v_dt[1]); 
            }

        }
          // compute mean and std 
        pair<double, double> h_trans = getMeanStd(hybrid_et); 
        pair<double, double> h_rot = getMeanStd(hybrid_er); 
        pair<double, double> t_trans = getMeanStd(epnp_et); 
        pair<double, double> t_rot = getMeanStd(epnp_er); 

        cout<<"nosie: "<<v_noise[i]<<" epnp mean_re: "<<t_rot.first<<" hybrid mean_re: "<<h_rot.first<<endl
            << "epnp mean_te: "<<t_trans.first<<" hybrid mean_te: "<<h_trans.first<<endl; 

        ouf_epnp<<v_noise[i]<<" \t "<<t_trans.first<<" \t "<<t_trans.second<<" \t "<<t_rot.first<<" \t "<<t_rot.second<<endl; 
        ouf_hybrid<<v_noise[i]<<" \t "<<h_trans.first<<" \t "<<h_trans.second<<" \t "<<h_rot.first<<" \t "<<h_rot.second<<endl; 
    }
    ouf_epnp.close(); 
    ouf_hybrid.close(); 
}

bool run_once(double noise, int num, vector<double>& v_dR, vector<double>& v_dt)
{
    SimCorr sim; 

    //initialize random seed
    initializeRandomSeed();

    //set experiment parameters
    // double noise = 0.; //1.; // 0.0;
    double outlierFraction = 0.0;
    size_t numberPoints = 1000;

    //create a random viewpoint pose
    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    Eigen::Matrix3d Rij = rotation;
    Vector3d tij = position; 
    Matrix3d Rji = Rij.transpose(); 
    Vector3d tji = -Rij.transpose()*tij; 

    // cout <<"gt rotation: "<<endl<<rotation<<endl; 

    //create a fake central camera
    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem( camOffsets, camRotations );

    //derive correspondences based on random point-cloud
    bearingVectors_t bearingVectors;
    points_t points;
    std::vector<int> camCorrespondences; //unused in the central case!
    Eigen::MatrixXd gt(3,numberPoints);
    generateRandom2D3DCorrespondences(
      position, rotation, camOffsets, camRotations, numberPoints, /*noise*/ 0, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

    vector<pair<Vector3d, Vector3d>> corrs; 
    for(int i=0; i<points.size(); i++){
        Vector3d pti = points[i]; 
        Vector3d ptj = bearingVectors[i]; 
        double pz = pti.z(); 
        pti = pti/pz;
        pti(2) = pz; 
        ptj = ptj/ptj(2); 
        corrs.push_back(make_pair(pti, ptj));  
    }

    // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 
    Matrix3d Rij_e, Rji_e; //, Rij_e_t;
    Vector3d tij_e, tji_e; // tij_e_t; 

    int cnt_3d = num; 
    OptSolver opt_solver; 

    {
        // cout<<"times: "<<i<<endl; 
        vector<pair<Vector3d, Vector3d>> corrs_noise = sim.addNoise3D2D(corrs, noise/460); 

        // vector<pair<Vector3d, Vector3d>> in_2d = getN(corrs_noise, cnt_2d); 
        vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs_noise, cnt_3d); 

        // 3d-2d result 
        me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
        // print_err("opencv epnp: ", Rij_e, tij_e, Rij, tij); 
        
        Matrix3d dR = Rij.transpose()*Rij_e;
        if(computeAngle(dR) > 1.){
            cout<<"what? computeAngle(dR): "<<computeAngle(dR)<<endl; 
            return false; 
        }

        Vector3d dt = tij_e - tij; 

        v_dR.push_back(computeAngle(dR));  
        v_dt.push_back(dt.norm()); 

        //derive correspondences based on random point-cloud
        bearingVectors_t points;
        bearingVectors_t bearingVectors1; 
        bearingVectors_t bearingVectors2;

        for(int j=0; j<in_3d.size(); j++){
            Vector3d pi = in_3d[j].first;
            pi(0) = pi(0)*pi(2); 
            pi(1) = pi(1)*pi(2);  
            Vector3d pj = in_3d[j].second; 
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
        st.solveTCeres(in_3d, Rij_e, tij_e); 
        dt = tij_e - tij; 

        v_dR.push_back(computeAngle(dR));  
        v_dt.push_back(dt.norm()); 

    }
    return true; 
}

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" rotation err: "<<computeAngle(dR)<<endl; 
}