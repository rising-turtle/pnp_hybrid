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
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>
#include "random_generators.hpp"
#include "experiment_helpers.hpp"

using namespace Eigen; 
using namespace std; 
using namespace opengv; 

void compare_hybrid(); 
void test_opengv_abs();
void test_opengv_rel();  
void print_err(string pre, Matrix<double, 3, 4>& Transformation, Matrix3d& Rg, Vector3d& tg);
void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 
void test_opt(); 

vector<pair<Vector3d, Vector3d>> combine(bearingVectors_t& pts, Eigen::MatrixXd& gt){
    vector<pair<Vector3d, Vector3d>> ret; 

    for(int col = 0; col < pts.size(); col++){
        Vector3d pti = gt.col(col); 
        Vector3d ptj = pts[col]; 
        pti(0) = pti(0)/pti(2); 
        pti(1) = pti(1)/pti(2); 
        ptj = ptj/ptj.z(); 
        ret.push_back(make_pair(pti, ptj)); 
    }
    return ret; 
}


int main(int argc, char* argv[])
{
    // test_opengv_rel();
    compare_hybrid(); 
    // test_opengv_abs(); 

    // test_opt(); 

    return 0; 
}

void test_opt(){
    ifstream inf("out.log");
    int N = 10;
    string line;
    double x1,y1,z1, x2,y2,z2; 
    vector<pair<Vector3d, Vector3d>> in_3d ; 
    while(inf.good()){
        stringstream ss;  
        getline(inf, line); 
        if(line.empty()) break; 
        ss << line; 
        ss >> x1 >> y1 >> z1 >> x2 >> y2 >> z2; 
        in_3d.push_back(make_pair(Vector3d{x1, y1, z1}, Vector3d{x2, y2, z2})); 
        cout << x1 << " "<< y1<<" "<< z1<<" "<<x2<<" "<<y2<<" "<<z2<<endl;
    }
    Matrix3d Rij_e; 
    Rij_e<<  0.87697,  0.311094, -0.366257,
            -0.211306,  0.934174,  0.287522,
            0.431594, -0.174756,  0.884979;

    // Rij_e <<  0.930365,   0.350576,   0.107321,
    //         -0.360711 ,  0.927647,   0.096738,
    //         -0.0656419,  -0.128713,   0.989507; 
    Vector3d tij_e{-0.528285, -1.08058, -0.619827};
    // Vector3d tij_e{-0.143357, -0.144381,  -1.12334};

    OptSolver opt_solver; 
    opt_solver.solveCeres(in_3d, Rij_e, tij_e); 

}

void compare_hybrid()
{
    SimCorr sim; 

    // stringstream ss; 
    // double roll(3.); // 7 
    // double yaw(1.); // 10
    // double pitch(2.);  // 5
    // cout<<" ground truth yaw: "<<yaw<<" pitch: "<<pitch<<" roll: "<<roll<<endl;

    // Eigen::AngleAxisd rollAngle(D2R(roll), Eigen::Vector3d::UnitX());
    // Eigen::AngleAxisd yawAngle(D2R(yaw), Eigen::Vector3d::UnitZ());
    // Eigen::AngleAxisd pitchAngle(D2R(pitch), Eigen::Vector3d::UnitY());
    // // Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
    // Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle; 
    // Eigen::Matrix3d Rij = q.matrix();

    // // cout<<"ground truth Rij: "<<endl<<Rij<<endl; 
    
    // Vector3d tij(0.2, 0.05, 0.3); 

    // Matrix3d Rji = Rij.transpose(); 
    // Vector3d tji = -Rij.transpose()*tij; 

    //initialize random seed
    initializeRandomSeed();

    //set experiment parameters
    double noise = 0.; //1.; // 0.0;
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
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
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

    // vector<pair<Vector3d, Vector3d>> corrs = sim.find_corrs(Rji, tji);
    
    // for(int k=0; k<corrs.size(); k++){
    //     Vector3d pi = corrs[k].first;
    //     pi(0) = pi(0)*pi(2); 
    //     pi(1) = pi(1)*pi(2);  
    //     Vector3d pj = corrs[k].second;
    //     Vector3d pj_hat = Rji*pi + tji; 
    //     pj_hat = pj_hat/pj_hat.z(); 
    //     cout <<"pi: "<<pi.transpose()<<"pj: "<<pj.transpose()<<" pj_hat: "<<pj_hat.transpose()<<endl;  
    //     if(k >= 7 ) break; 
    // }

    // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 
    Matrix3d Rij_e, Rji_e; //, Rij_e_t;
    Vector3d tij_e, tji_e; // tij_e_t; 

    Matrix3d R2; //, Rij_e_t;
    Vector3d t2; // tij_e_t; 

    int cnt_2d = 40; 
    int cnt_3d = 10; 
    OptSolver opt_solver; 

    // 3d-2d 
    cv::Mat rvec, tvec; 
    Matrix3d tmpR; 
    Vector3d tmpt; 

    for(int i=0; i<1; i++){

        // cout<<"times: "<<i<<endl; 
        vector<pair<Vector3d, Vector3d>> corrs_noise = sim.addNoise3D2D(corrs); 

        // vector<pair<Vector3d, Vector3d>> in_2d = getN(corrs_noise, cnt_2d); 
        vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs_noise, cnt_3d); 

        // 3d-2d result 
        me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
        print_err("opencv epnp: ", Rij_e, tij_e, Rij, tij); 
        
        // run 2D-2D result, eigen solver 
        // me.solvePNP_2D_2D(in_3d, tmpR, tmpt); 
        
        // Matrix3d dR1 = Rij.transpose()*tmpR;
        // cout<<"eight_pt err: "<<computeAngle(dR1)<<endl; 

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
        
        // Eigen::Matrix3d rotation = Eigen::Matrix3d::Identity(); 
        //create a central absolute adapter
        // absolute_pose::CentralAbsoluteAdapter adapter_abs(
        //     bearingVectors2,
        //     points,
        //     rotation);

        // transformation_t epnp_transformation = absolute_pose::epnp(adapter_abs);
        // print_err("opengv epnp : ", epnp_transformation, Rij, tij);
        // print_err("opengv epnp : ", epnp_transformation, Rji, tji);

        relative_pose::CentralRelativeAdapter adapter_rbs(
              bearingVectors1,
              bearingVectors2,
              Rij_e ); // rotation);

        // first use eight pts to compute initial rotation 

        Rij_e =  relative_pose::eigensolver(adapter_rbs);  
        // Rij_e = Rji_e.transpose(); 
        
        Matrix3d dR = Rij.transpose()*Rij_e;
        cout<<"eigensolver_rotation err: "<<computeAngle(dR)<<endl; 
        
        st.solveTCeres(in_3d, Rij_e, tij_e); 
        // optimize the [Rij_e] and [tij_e]
        print_err("initial hybrid: ", Rij_e, tij_e, Rij, tij);

        ofstream ouf("out.log"); 
        for(int i=0; i<in_3d.size(); i++){
            ouf<<in_3d[i].first.transpose()<<" "<<in_3d[i].second.transpose()<<endl; 
        }
        cout<<"Rij_e: "<<endl<<Rij_e<<endl<<" tij: "<<tij.transpose()<<endl; 
        opt_solver.solveCeres(in_3d, Rij_e, tij_e); 
        print_err("hybrid after opt: ", Rij_e, tij_e, Rij, tij);

        cout<<endl<<endl;
    }
    return ; 
}

void test_opengv_rel()
{
    // initialize random seed
    initializeRandomSeed();

    //set experiment parameters
    double noise = 0.0;
    double outlierFraction = 0.0;
    size_t numberPoints = 10;

    //generate a random pose for viewpoint 1
    translation_t position1 = Eigen::Vector3d::Zero();
    rotation_t rotation1 = Eigen::Matrix3d::Identity();


    double roll(7.); // 7 
    double yaw(10.); // 10
    double pitch(5.);  // 5
    cout<<" ground truth yaw: "<<yaw<<" pitch: "<<pitch<<" roll: "<<roll<<endl;

    Eigen::AngleAxisd rollAngle(D2R(roll), Eigen::Vector3d::UnitX());
    Eigen::AngleAxisd yawAngle(D2R(yaw), Eigen::Vector3d::UnitZ());
    Eigen::AngleAxisd pitchAngle(D2R(pitch), Eigen::Vector3d::UnitY());
    // Eigen::Quaternion<double> q = rollAngle * yawAngle * pitchAngle;
    Eigen::Quaternion<double> q = yawAngle * pitchAngle * rollAngle; 
    Eigen::Matrix3d Rij = q.matrix();

    // cout<<"ground truth Rij: "<<endl<<Rij<<endl; 
    
    Vector3d tij(0.2, 0.05, 0.3); 

    Matrix3d Rji = Rij.transpose(); 
    Vector3d tji = -Rij.transpose()*tij; 

    //generate a random pose for viewpoint 2
    translation_t position2 = generateRandomTranslation(2.0);
    rotation_t rotation2 = generateRandomRotation(0.5);

    //create a fake central camera
    translations_t camOffsets;
    rotations_t camRotations;
    generateCentralCameraSystem( camOffsets, camRotations );

    //derive correspondences based on random point-cloud
    bearingVectors_t bearingVectors1;
    bearingVectors_t bearingVectors2;
    std::vector<int> camCorrespondences1; //unused in the central case
    std::vector<int> camCorrespondences2; //unused in the central case
    Eigen::MatrixXd gt(3,numberPoints);
    generateRandom2D2DCorrespondences(
      position1, rotation1, position2, rotation2,
      camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors1, bearingVectors2,
      camCorrespondences1, camCorrespondences2, gt );

    //Extract the relative pose
    translation_t position; rotation_t rotation;
    extractRelativePose(
      position1, position2, rotation1, rotation2, position, rotation );


    rotation_t tmp = Matrix3d::Identity(); 
    // translation_t tmp_t = Vector3d::Identity(); 

    //create a central relative adapter
    relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2,
      tmp);

    rotation_t eigensolver_rotation = relative_pose::eigensolver(adapter);

    Matrix3d dR = rotation.transpose()*eigensolver_rotation;

    cout<<"eigensolver_rotation err: "<<computeAngle(dR)<<endl; 

}


void test_opengv_abs()
{
    /*
    //initialize random seed
    initializeRandomSeed();

    //set experiment parameters
    double noise = 1.; //1.; // 0.0;
    double outlierFraction = 0.0;
    size_t numberPoints = 100;

    //create a random viewpoint pose
    translation_t position = generateRandomTranslation(2.0);
    rotation_t rotation = generateRandomRotation(0.5);

    cout <<"gt rotation: "<<endl<<rotation<<endl; 

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
      position, rotation, camOffsets, camRotations, numberPoints, noise, outlierFraction,
      bearingVectors, points, camCorrespondences, gt );

    
    rotation_t tmp = Matrix3d::Identity(); 

    //create a central absolute adapter
    absolute_pose::CentralAbsoluteAdapter adapter(
      bearingVectors,
      points,
      tmp );

    transformation_t epnp_transformation = absolute_pose::epnp(adapter);
    print_err("epnp : ", epnp_transformation, rotation, position);
    cout <<"opengv epnp: "<<endl<<epnp_transformation.block<3,3>(0,0)<<endl; 

    std::vector<int> indices6 = getNindices(6);
    transformation_t epnp_transformation_6 =
      absolute_pose::epnp( adapter, indices6 );
    print_err("epnp_6 : ", epnp_transformation_6, rotation, position);    

    transformation_t nonlinear_transformation = absolute_pose::optimize_nonlinear(adapter);
    print_err("nonlinear_transformation : ", nonlinear_transformation, rotation, position);
  
    vector<pair<Vector3d, Vector3d>> in_3d = combine(bearingVectors, gt); 
    MotionEstimator me;
    Matrix3d Rij_e; //, Rij_e_t;
    Vector3d tij_e; // tij_e_t; 
    // 3d-2d result 
    me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
    cout <<"Rij_e: "<<endl<<Rij_e<<endl; 
    print_err("EPNP: ", Rij_e, tij_e, rotation, position); */
    return ; 
}





void print_err(string pre, Matrix<double, 3, 4>& Transformation, Matrix3d& Rg, Vector3d& tg)
{

    Matrix3d Re = Transformation.block<3,3>(0, 0);
    Vector3d te = Transformation.block<3,1>(0, 3);
    Vector3d dt = tg - te; 
    Matrix3d dR = Rg.transpose()*Re;

    cout<<pre<<" trans err: "<<dt.norm()<<" ros err: "<<computeAngle(dR)<<endl; 
}

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
    Vector3d dt = tg - te; 
    Matrix3d dR = Rg.transpose()*Re;

    cout<<pre<<" trans err: "<<dt.norm()<<" rotation err: "<<computeAngle(dR)<<endl; 
}