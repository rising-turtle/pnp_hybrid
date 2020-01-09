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
#include <sstream>
#include <fstream>

using namespace Eigen; 
using namespace std; 

// void test_hybrid(SimCorr& ); 
void test_2d_2d( vector<pair<Vector3d, Vector3d>>& corrs_noise, Matrix3d& Rij_e, Vector3d& tij_e, cv::Mat& mask); 
void test_3d_2d( vector<pair<Vector3d, Vector3d>>& corrs_noise, Matrix3d& Rij_e, Vector3d& tij_e, cv::Mat& mask);

pair<double, double> run_once_3d_2d(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt);

pair<double, double> run_once_hybrid(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt);

vector<double> run_once_together(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt);

void run_monte_carlo(vector<int> v_cnt_3d, int cnt_2d = 40, int TIMES = 10); 

int main(int argc, char* argv[])
{
    // vector<int> v3d{4, 7, 10, 15, 20, 25, 30, 35, 40}; 
    vector<int> v3d{30}; 
    run_monte_carlo(v3d, 40, 50); 

    return 0; 
}

void run_monte_carlo(vector<int> v_cnt_3d, int cnt_2d, int TIMES)
{
    SimCorr sim; 

    stringstream ss; 
    ss<<"output_cnt_2d_"<<cnt_2d<<".log"; 
    ofstream ouf(ss.str().c_str()); 

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

    ouf<<"Num_of_3d "<<"\t"<<" hybrid_trans "<<"\t"<<" hybrid_rot "<<"\t"<<" 3d_2d_trans "<<"\t"<<" 3d_2d_rot "<<
    "\t"<<" hybrid_trans "<<"\t"<<" hybrid_rot "<<"\t"<<" 3d_2d_trans "<<"\t"<<" 3d_2d_rot"<<endl; 

    for(int j=0; j<v_cnt_3d.size(); j++){
        int cnt_3d = v_cnt_3d[j]; 
        vector<double> h_et, h_ea, t_et, t_ea; 
        for(int cnt = 0; cnt < TIMES; cnt++){
            vector<pair<Vector3d, Vector3d>> corrs_noise = sim.add_noise(corrs); 
            // pair<double, double> hybrid_err = run_once_hybrid(corrs_noise, cnt_3d, cnt_2d, Rij, tij); 
            // pair<double, double> t32_err = run_once_3d_2d(corrs_noise, cnt_3d, cnt_2d, Rij, tij);
            // h_et.push_back(hybrid_err.first); h_ea.push_back(hybrid_err.second); 
            // t_et.push_back(t32_err.first); t_ea.push_back(t32_err.second); 

            vector<double> err = run_once_together(corrs_noise, cnt_3d, cnt_2d, Rij, tij); 
            h_et.push_back(err[0]); h_ea.push_back(err[1]); 
            if(err[2] < 2.)
                t_et.push_back(err[2]); 
            if(err[3] < D2R(30))
                t_ea.push_back(err[3]); 

            // cout<<"cnt: "<<cnt<<" hybrid error: trans: "<<hybrid_err.first<<" angle: "<<hybrid_err.second<<" 3d-2d error: trans: "<<
            //     t32_err.first<<" angle: "<< t32_err.second<<endl;

            cout<<"cnt: "<<cnt<<" hybrid error: trans: "<<err[0]<<" angle: "<<err[1]<<" 3d-2d error: trans: "<<
                 err[2]<<" angle: "<< err[3]<<endl;
        }

        // compute mean and std 
        pair<double, double> h_trans = getMeanStd(h_et); 
        pair<double, double> h_rot = getMeanStd(h_ea); 
        pair<double, double> t_trans = getMeanStd(t_et); 
        pair<double, double> t_rot = getMeanStd(t_ea); 

        ouf<<cnt_3d<<" \t "<<h_trans.first<<" \t "<<h_rot.first<<" \t "<<t_trans.first<<" \t "<<t_rot.first<<" \t "<<
            h_trans.second<<" \t "<<h_rot.second<<" \t "<<t_trans.second<<" \t "<<t_rot.second<<endl;
       //  vector<pair<double, double>> ret(4); 
        // ret[0] = h_trans; ret[1] = h_rot; ret[2] = t_trans; ret[3] = t_rot; 
    }
    ouf.close(); 
    return ; 
}

vector<double> run_once_together(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt)
{
    Matrix3d Rij_e_h, Rij_e_t;
    Vector3d tij_e_h, tij_e_t; 
    cv::Mat mask, mask_t; 
    test_2d_2d(corres, Rij_e_h, tij_e_h, mask); 

    // with the speficied number of features 
    vector<pair<Vector3d, Vector3d>> inliers = getInliers(corres, mask); 

    // from the 2d inliers find out 3d inliers 
    test_3d_2d(inliers, Rij_e_t, tij_e_t, mask_t); 

    vector<pair<Vector3d, Vector3d>> inliers_3d = getInliersIndex(inliers, mask_t); 
    vector<pair<Vector3d, Vector3d>> in_2d = getN(inliers, cnt_2d); 
    vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers_3d, cnt_3d); 

    // update Rij_e using the specified number of 2d features 
    MotionEstimator me;
    me.solvePNP_2D_2D(in_2d, Rij_e_h, tij_e_h);    

    // compute 3d-2d transformation 
    me.solvePNP_3D_2D(in_3d, Rij_e_t, tij_e_t); 

    // estimate translation 
    SolveTranslate st; 
    Eigen::Vector3d ntij = tij_e_h; 
    tij_e_h = tij_e_t; // use 3d-2d results as initial value 
    // st.solveTCeres(in_3d, Rij_e_h, tij_e_h); 
    // st.solveTProjCeres(in_3d, Rij_e_h, tij_e_h);
    // st.solveTScaleCeres(in_3d, Rij_e_h, ntij, tij_e_h); 
    // st.solveTProjScaleCeres(in_3d, Rij_e_h, ntij, tij_e_h); 
    // st.solveTCeresWithPt(in_3d, Rij_e_h, tij_e_h); 
    st.solveTProjCeresWithPt(in_3d, Rij_e_h, tij_e_h); 

    // compute error 
    Matrix3d dR_h = Rij_gt.transpose()*Rij_e_h; 
    Vector3d dt_h = tij_gt - tij_e_h; 

    double et_h = dt_h.norm(); 
    double ea_h = computeAngle(dR_h); 

 

    Matrix3d dR_t = Rij_gt.transpose()*Rij_e_t; 
    Vector3d dt_t = tij_gt - tij_e_t; 

    double et_t = dt_t.norm(); 
    double ea_t = computeAngle(dR_t); 


    cout<<"T_3d_2d err: "<<sum_error(in_3d, Rij_e_h, tij_e_t)<<" T_gt err: "<<sum_error(in_3d, Rij_e_h, tij_gt)<<
        " T_hybrid err: "<<sum_error(in_3d, Rij_e_h, tij_e_h)<<endl; 


    vector<double> ret{et_h, ea_h, et_t, ea_t}; 

    return ret; 
}


pair<double, double> run_once_hybrid(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt)
{

    Matrix3d Rij_e;
    Vector3d tij_e; 
    cv::Mat mask; 
    test_2d_2d(corres, Rij_e, tij_e, mask); 

    // with the speficied number of features 
    vector<pair<Vector3d, Vector3d>> inliers = getInliers(corres, mask); 
    vector<pair<Vector3d, Vector3d>> in_2d = getN(inliers, cnt_2d); 
    vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers, cnt_3d); 

    // update Rij_e using the specified number of 2d features 
    mask = cv::Mat();
    // test_2d_2d(in_2d, Rij_e, tij_e, mask); 
    MotionEstimator me;
    me.solvePNP_2D_2D(in_2d, Rij_e, tij_e);    

    // estimate translation 
    SolveTranslate st; 
    st.solveTCeres(in_3d, Rij_e, tij_e); 

    // compute error 
    Matrix3d dR = Rij_gt.transpose()*Rij_e; 
    Vector3d dt = tij_gt - tij_e; 

    double et = dt.norm(); 
    double ea = computeAngle(dR); 
    return make_pair(et, ea); 
}



pair<double, double> run_once_3d_2d(vector<pair<Vector3d, Vector3d>>& corres, int cnt_3d, int cnt_2d, 
    Matrix3d& Rij_gt, Vector3d& tij_gt)
{
    Matrix3d Rij_e;
    Vector3d tij_e; 
    cv::Mat mask;

    test_3d_2d(corres, Rij_e, tij_e, mask); 

    // with the speficied number of features 
    cout<<"corres.size(): "<<corres.size()<<" mask.rows: "<<mask.rows<<endl;
    vector<pair<Vector3d, Vector3d>> inliers = getInliersIndex(corres, mask); 
    vector<pair<Vector3d, Vector3d>> in_3d = getN(inliers, cnt_3d); 

    mask = cv::Mat(); 
    // test_3d_2d(in_3d, Rij_e, tij_e, mask); 
    
    MotionEstimator me;
    me.solvePNP_3D_2D(in_3d, Rij_e, tij_e);    

    // compute error 
    Matrix3d dR = Rij_gt.transpose()*Rij_e; 
    Vector3d dt = tij_gt - tij_e; 

    double et = dt.norm(); 
    double ea = computeAngle(dR); 
    return make_pair(et, ea); 
}

void test_2d_2d( vector<pair<Vector3d, Vector3d>>& corrs_noise, Matrix3d& Rij_e, Vector3d& tij_e, cv::Mat& mask)
{
    MotionEstimator me;
    // 2d-2d  
    me.solveRelativeRT(corrs_noise, Rij_e, tij_e, mask); 
    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 
    // cout << "to Euler angles:" << endl;
    // cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    // cout <<"main.cpp: 2d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 
}


void test_3d_2d(vector<pair<Vector3d, Vector3d>>& corrs_noise, Matrix3d& Rij_e, Vector3d& tij_e, cv::Mat& mask)
{
    MotionEstimator me;

    // 3d-2d 
    me.solveRelativeRT_PNP(corrs_noise, Rij_e, tij_e, mask);
    Vector3d ea = Rij_e.eulerAngles(2, 1, 0); 

    // cout << "to Euler angles:" << endl;
    // cout << "yaw: "<<R2D(ea.x()) <<" pitch: "<<R2D(ea.y())<<" roll: "<<R2D(ea.z()) << endl;
    // cout <<"main.cpp: 3d-2d estimate Rij: "<<endl<<Rij_e<<endl<<"tij: "<<tij_e.transpose()<<endl; 
}
