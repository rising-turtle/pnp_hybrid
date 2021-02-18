/*
    Feb. 9th, 2021, He Zhang, fuyinzh@gmail.com 

    monte carlo against different number of points   

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

void print_err(string pre, Matrix<double, 3, 4>& Transformation, Matrix3d& Rg, Vector3d& tg);
void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 

bool run_once(double noise, int cnt_3d, int cnt_2d, vector<double>& dR, vector<double>& dt, bool use_optimization=false); 

void run_monte_carlo( vector<int> v_cnt_3d, double noise = 1., int cnt_2d = 30,  int TIMES= 1000) ;

int main(int argc, char* argv[])
{
    // set up camera intrinsic 
    CX = 960; CY = 650; 
    FX = FY = 1460; 

    vector<int> v_cnt_3d{ 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    // vector<int> v_cnt_3d{ 10, 20, 28};
    run_monte_carlo(v_cnt_3d, 2.0, 30, 2000);
    // vector<double> v_dR, v_dt; 
    // run_once(1., 10, 30, v_dR, v_dt); 
    return 0; 
}


void run_monte_carlo( vector<int> v_cnt_3d, double noise, int cnt_2d,  int TIMES)
{
    // ofstream ouf_epnp("PT_EPNP_NOISE_20.log"); 
    // ofstream ouf_hybrid("PT_HYBRID_NOISE_20.log"); 

    ofstream ouf("PT_NOISE_20.log"); 
    ouf <<"number of 3D poitns: [mean_t] [std_t] [mean_r] [std_r]"<<endl; 
    vector<string> methods{"EPNP", "HYBRID", "HYBRID-GN"}; // , "UPNP"

    // ofstream ouf_epnp("tmp_epnp.log"); 
    // ofstream ouf_hybrid("tmp_hybrid.log"); 

    int N = methods.size(); 
    for(int n = 0; n<N ; n++)
        ouf << methods[n]<<" ";
    ouf<<std::fixed<<endl; 

    for(int i=0; i<v_cnt_3d.size(); i++){

        // vector<double> hybrid_et, hybrid_er, epnp_et, epnp_er; 
        vector<vector<double> > vv_er(N);
        vector<vector<double> > vv_et(N);

        for(int k=0; k < TIMES; ){
            vector<double> v_dr, v_dt; 
            if(run_once(noise, v_cnt_3d[i], cnt_2d, v_dr, v_dt, true)){
                k++; 
                for(int n=0; n<N; n++){
                    vv_er[n].push_back(v_dr[n]);
                    vv_et[n].push_back(v_dt[n]); 
                }
                // epnp_er.push_back(v_dr[0]); 
                // epnp_et.push_back(v_dt[0]);
                // hybrid_er.push_back(v_dr[1]);
                // hybrid_et.push_back(v_dt[1]); 
            }
        }

        vector<double> v_mean_r(N); 
        vector<double> v_std_r(N); 
        vector<double> v_mean_t(N); 
        vector<double> v_std_t(N);
        for(int n=0; n<N; n++){
            pair<double, double> trans = getMeanStd(vv_et[n]); 
            v_mean_t[n] = trans.first; 
            v_std_t[n] = trans.second; 
            pair<double, double> rot = getMeanStd(vv_er[n]);
            v_mean_r[n] = rot.first; 
            v_std_r[n] = rot.second;  
        } 
          // compute mean and std 
        // pair<double, double> h_trans = getMeanStd(hybrid_et); 
        // pair<double, double> h_rot = getMeanStd(hybrid_er); 
        // pair<double, double> t_trans = getMeanStd(epnp_et); 
        // pair<double, double> t_rot = getMeanStd(epnp_er); 


        ouf<<v_cnt_3d[i]<<" ";
        for(int n=0; n<N; n++){
            cout<<"cnt_3d: "<<v_cnt_3d[i]<<" "<<methods[n]<<" mean_te: "<<v_mean_t[n]<<" mean_re: "<<v_mean_r[n]<<endl; 
            ouf<<v_mean_t[n]<<" "<<v_std_t[n]<<" "<<v_mean_r[n]<<" "<<v_std_r[n]<<" ";
        }
        ouf<<endl;

        // ouf_epnp<<v_cnt_3d[i]<<" \t "<<t_trans.first<<" \t "<<t_trans.second<<" \t "<<t_rot.first<<" \t "<<t_rot.second<<endl; 
        // ouf_hybrid<<v_cnt_3d[i]<<" \t "<<h_trans.first<<" \t "<<h_trans.second<<" \t "<<h_rot.first<<" \t "<<h_rot.second<<endl; 
    }
    // ouf_epnp.close(); 
    // ouf_hybrid.close(); 
    ouf.close(); 
}


bool run_once(double noise, int cnt_3d, int cnt_2d, vector<double>& v_dR, vector<double>& v_dt, bool use_optimization)
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

    OptSolver opt_solver; 

    {
        // cout<<"times: "<<i<<endl; 
        vector<pair<Vector3d, Vector3d>> corrs_noise = sim.addNoise3D2D(corrs, noise/FX); 

        vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs_noise, cnt_3d); 

        vector<pair<Vector3d, Vector3d>> in_2d = in_3d;

        vector<pair<Vector3d, Vector3d>> in_2d_tmp;
        if(cnt_2d > cnt_3d){
            in_2d_tmp = getN(corrs_noise, cnt_2d-cnt_3d); 
            in_2d.insert(in_2d.end(), in_2d_tmp.begin(), in_2d_tmp.end()); 
        }

        // 3d-2d result 
        me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
        print_err("opencv epnp: ", Rij_e, tij_e, Rij, tij); 
        
        Matrix3d dR = Rij.transpose()*Rij_e;
        if(computeAngle(dR) > 1.){
            // cout<<"what? computeAngle(dR): "<<computeAngle(dR)<<endl; 
            return false; 
        }

        Vector3d dt = tij_e - tij; 
        if(dt.norm() > 1.){
            return false; 
        }

        v_dR.push_back(computeAngle(dR));  
        v_dt.push_back(dt.norm()); 

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
        st.solveTCeres(in_3d, Rij_e, tij_e); 
        dt = tij_e - tij; 

        cout<<"Hybrid PnP dR: "<<computeAngle(dR)<<" dt: "<<dt.norm()<<endl; 
        v_dR.push_back(computeAngle(dR));  
        v_dt.push_back(dt.norm()); 

        Matrix3d RR1 = Rij_e; Matrix3d RR2 = Rij_e; 
        Vector3d tt1 = tij_e; Vector3d tt2 = tij_e; 
        // opt_solver.solveCeres(in_3d, RR1, tt1); 
        // dR = Rij.transpose()*RR1; 
        // dt = tt1 - tij; 
        // cout<<"Opti_1 dR: "<<computeAngle(dR)<<" dt: "<<dt.norm()<<endl; 
        if(use_optimization){
            opt_solver.solveCeresHybrid(in_3d, in_2d_tmp, RR2, tt2); 
            dR = Rij.transpose()*RR2; 
            dt = tt2 - tij; 
            cout<<"Opti_2 dR: "<<computeAngle(dR)<<" dt: "<<dt.norm()<<endl;
            v_dR.push_back(computeAngle(dR));  
            v_dt.push_back(dt.norm());  
        }

        // // try UPNP
        // rotation_t tmp = Matrix3d::Identity(); 

        // //create a central absolute adapter
        // absolute_pose::CentralAbsoluteAdapter adapter_abs(
        //   bearingVectors2,
        //   points,
        //   tmp );

        // transformations_t upnp_transformations = absolute_pose::upnp(adapter_abs);
        // transformation_t upnp_transformation = upnp_transformations[0]; 
        // // for(int k=0; k < upnp_transformations.size(); k++)
        //     // print_err("upnp : ", upnp_transformations[k], Rij, tij);
        // // cout <<"opengv upnp: "<<endl<<upnp_transformation.block<3,3>(0,0)<<endl; 
        // Matrix3d R_upnp = upnp_transformation.block<3,3>(0, 0);
        // Vector3d t_upnp = upnp_transformation.block<3,1>(0, 3);
        // dt = tij - t_upnp; 
        // dR = Rij.transpose()*R_upnp;

        // v_dR.push_back(computeAngle(dR));  
        // v_dt.push_back(dt.norm()); 

    }
    return true; 
}

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" rotation err: "<<computeAngle(dR)<<endl; 
}

void print_err(string pre, Matrix<double, 3, 4>& Transformation, Matrix3d& Rg, Vector3d& tg)
{

    Matrix3d Re = Transformation.block<3,3>(0, 0);
    Vector3d te = Transformation.block<3,1>(0, 3);
    Vector3d dt = tg - te; 
    Matrix3d dR = Rg.transpose()*Re;

    cout<<pre<<" trans err: "<<dt.norm()<<" ros err: "<<computeAngle(dR)<<endl; 
}