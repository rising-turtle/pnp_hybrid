/*
    Jan. 6th, 2020, He Zhang, hzhang8@vcu.edu 

    entry file for simulation 

*/

#include "sim_corres.h"
#include "solve_5pts.h"
#include "solve_translate.h"
#include "utility.h"
#include <Eigen/Core>
#include <cmath>

using namespace Eigen; 
using namespace std; 

#define D2R(d) (((d)/180.)*M_PI)
#define R2D(r) (((r)/M_PI)*180.)


// void test_hybrid(SimCorr& ); 
void test_2d_2d( vector<pair<Vector3d, Vector3d>>& ); 
void test_3d_2d( vector<pair<Vector3d, Vector3d>>& );

int main(int argc, char* argv[])
{
    SimCorr sim; 

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

    SolveTranslate st; 
    vector<pair<Vector3d, Vector3d>> inliers = getN(corrs_noise, 20);

    st.solveTCeres(inliers, Rij, tij); 
    cout<<" main.cpp: hybrid 3d-2d tij: "<<tij.transpose()<<endl;

    return 0; 
}

void test_2d_2d( vector<pair<Vector3d, Vector3d>>& corrs_noise)
{
    MotionEstimator me;
    Matrix3d Rij_e;
    Vector3d tij_e; 
    cv::Mat mask; 
    SolveTranslate st; 

    // 2d-2d  
    me.solveRelativeRT(corrs_noise, Rij_e, tij_e, mask); 
    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 
    cout << "to Euler angles:" << endl;
    cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    cout <<"main.cpp: 2d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 
}


void test_3d_2d(vector<pair<Vector3d, Vector3d>>& corrs_noise )
{
    MotionEstimator me;
    Matrix3d Rij_e;
    Vector3d tij_e; 
    cv::Mat mask; 
    SolveTranslate st; 

    // 3d-2d 
    me.solveRelativeRT_PNP(corrs_noise, Rij_e, tij_e);
    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 

    cout << "to Euler angles:" << endl;
    cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    cout <<"main.cpp: 3d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 
}
