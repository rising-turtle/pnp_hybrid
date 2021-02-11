/*
    Feb. 10th, 2021, He Zhang, fuyinzh@gmail.com 

    monte carlo against different number of points with iphone's data    

*/

#include "solve_5pts.h"
#include "solve_translate.h"
#include "utility.h"
#include "opt_solver.h"
#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <opencv2/opencv.hpp>

#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/absolute_pose/methods.hpp>
#include <opengv/absolute_pose/CentralAbsoluteAdapter.hpp>

using namespace Eigen; 
using namespace std; 
using namespace opengv; 


string folder("/home/davidz/work/data/iphone_egomotion"); 
string move_folder_1("rot_yaw_4"); 
string move_folder_2("rot_yaw_7"); 

int COL = 1920; 
int ROW = 1440; 
double CX = 965.24;
double CY = 653.62; 
double FX = 1460.2;
double FY = 1460.2; 

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg); 

bool run_once(cv::Mat img_1, cv::Mat img_2, cv::Mat depth_img1, int cnt_3d, int cnt_2d, bool use_optimization=false); 

// void run_monte_carlo( vector<int> v_cnt_3d, int cnt_2d = 30,  int TIMES= 100) ;

void showMatchedFeatures(cv::Mat img_1, cv::Mat img_2, vector<pair<Vector3d, Vector3d>>& matches);

int main(int argc, char* argv[])
{
    // set up camera intrinsic 


    vector<int> v_cnt_3d{ 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30};
    // vector<int> v_cnt_3d{ 10, 20, 28};
    // run_monte_carlo(v_cnt_3d, 2.0, 30, 2000);
    // vector<double> v_dR, v_dt; 

    int file_num = 0;
    string fc_name_1 = folder+"/"+move_folder_1+"/color/"+to_string(200+file_num)+".png"; 
    string fc_name_2 = folder+"/"+move_folder_2+"/color/"+to_string(200+file_num)+".png";
    string fd_name_1 = folder+"/"+move_folder_1+"/depth/"+to_string(200+file_num)+".png";

    cv::Mat img_1 = cv::imread(fc_name_1.c_str(), -1); 
    cv::Mat img_2 = cv::imread(fc_name_2.c_str(), -1); 
    cv::Mat gray_1, gray_2; 
    cv::cvtColor(img_1, gray_1, CV_BGR2GRAY); 
    cv::cvtColor(img_2, gray_2, CV_BGR2GRAY); 

    cv::Mat dpt_1 = cv::imread(fd_name_1.c_str(), -1); 

    run_once(gray_1, gray_2, dpt_1, 10, 40); 
    return 0; 
}


// void run_monte_carlo( vector<int> v_cnt_3d, double noise, int cnt_2d,  int TIMES)
// {
//     ofstream ouf_epnp("PT_EPNP_NOISE_20.log"); 
//     ofstream ouf_hybrid("PT_HYBRID_NOISE_20.log"); 

//     // ofstream ouf_epnp("tmp_epnp.log"); 
//     // ofstream ouf_hybrid("tmp_hybrid.log"); 

//     for(int i=0; i<v_cnt_3d.size(); i++){

//         vector<double> hybrid_et, hybrid_er, epnp_et, epnp_er; 

//         for(int k=0; k < TIMES; ){
//             vector<double> v_dr, v_dt; 
//             if(run_once(noise, v_cnt_3d[i], cnt_2d, v_dr, v_dt)){
//                 k++; 
//                 epnp_er.push_back(v_dr[0]); 
//                 epnp_et.push_back(v_dt[0]);
//                 hybrid_er.push_back(v_dr[1]);
//                 hybrid_et.push_back(v_dt[1]); 
//             }

//         }
//           // compute mean and std 
//         pair<double, double> h_trans = getMeanStd(hybrid_et); 
//         pair<double, double> h_rot = getMeanStd(hybrid_er); 
//         pair<double, double> t_trans = getMeanStd(epnp_et); 
//         pair<double, double> t_rot = getMeanStd(epnp_er); 

//         cout<<"cnt_3d: "<<v_cnt_3d[i]<<" epnp mean_re: "<<t_rot.first<<" hybrid mean_re: "<<h_rot.first<<endl
//             << "epnp mean_te: "<<t_trans.first<<" hybrid mean_te: "<<h_trans.first<<endl; 

//         ouf_epnp<<v_cnt_3d[i]<<" \t "<<t_trans.first<<" \t "<<t_trans.second<<" \t "<<t_rot.first<<" \t "<<t_rot.second<<endl; 
//         ouf_hybrid<<v_cnt_3d[i]<<" \t "<<h_trans.first<<" \t "<<h_trans.second<<" \t "<<h_rot.first<<" \t "<<h_rot.second<<endl; 
//     }
//     ouf_epnp.close(); 
//     ouf_hybrid.close(); 
// }


bool inBorder(const cv::Point2f &pt)
{
    const int BORDER_SIZE = 7;
    int img_x = cvRound(pt.x);
    int img_y = cvRound(pt.y);
    return BORDER_SIZE <= img_x && img_x < COL - BORDER_SIZE && BORDER_SIZE <= img_y && img_y < ROW - BORDER_SIZE;
}


void reduceVector(vector<cv::Point2f> &v, vector<uchar> status)
{
    int j = 0;
    for (int i = 0; i < int(v.size()); i++)
        if (status[i])
            v[j++] = v[i];
    v.resize(j);
}


void rejectWithF( vector<cv::Point2f>& pts_1, vector<cv::Point2f>& pts_2)
{
    if (pts_1.size() >= 8)
    {
        vector<uchar> status;
        double F_THRESHOLD = 0.5;
        cv::findFundamentalMat(pts_1, pts_2, cv::FM_RANSAC, F_THRESHOLD, 0.99, status);
        reduceVector(pts_1, status);
        reduceVector(pts_2, status);
    }
}


void find3D2DMatch(const vector<cv::Point2f>& pts_1, const vector<cv::Point2f>& pts_2, 
    cv::Mat depth, vector<pair<Vector3d, Vector3d>>& in_3d, vector<pair<Vector3d, Vector3d>>& in_2d)
{
    double ui, vi, uj, vj; 
    double dpt_scale = 7.5; 
    double xi, yi, xj, yj; 
    in_3d.clear(); 
    in_2d.clear(); 
    for(int i=0; i < pts_1.size(); i++){

        ui = pts_1[i].x; 
        vi = pts_1[i].y; 

        uj = pts_2[i].x; 
        vj = pts_2[i].y; 

        xi = (ui - CX)/FX; 
        yi = (vi - CY)/FY; 
        xj = (uj - CX)/FX;
        yj = (vj - CY)/FY; 

        int ui_dpt = round(ui/dpt_scale); 
        int vi_dpt = round(vi/dpt_scale); 

        double dist_ui = ui_dpt - ui/dpt_scale; 
        double dist_vi = vi_dpt - vi/dpt_scale; 

        double d = depth.at<unsigned short>(vi_dpt, ui_dpt) * 0.001; 

        bool close = false; 
        if(fabs(dist_ui) <0.4 && fabs(dist_vi) < 0.4)
            close = true; 

        if(d >= 0.2 && d <= 4.5 && close){
            // valid 3d point 
            in_3d.push_back(make_pair(Vector3d{xi, yi, d}, Vector3d{xj, yj, 1.})); 
        }else{
            in_2d.push_back(make_pair(Vector3d{xi, yi, 1.}, Vector3d{xj, yj, 1.})); 
        }
    }
    return ; 
}

// bool run_once(double noise, int cnt_3d, int cnt_2d, vector<double>& v_dR, vector<double>& v_dt, bool use_optimization)
bool run_once(cv::Mat img_1, cv::Mat img_2, cv::Mat depth_img1, int cnt_3d, int cnt_2d, bool use_optimization)
{
    size_t numberPoints = 1000;
    vector<cv::Point2f> pts_1;
    double MIN_DIST = 20.;
    cv::goodFeaturesToTrack(img_1, pts_1, numberPoints, 0.01, MIN_DIST);

    vector<uchar> status;
    vector<float> err;
    vector<cv::Point2f> pts_2;
    cv::calcOpticalFlowPyrLK(img_1, img_2, pts_1, pts_2, status, err, cv::Size(21, 21), 3);

    for (int i = 0; i < int(pts_2.size()); i++)
        if (status[i] && !inBorder(pts_2[i]))
                status[i] = 0;

    reduceVector(pts_1, status); 
    reduceVector(pts_2, status); 
    rejectWithF(pts_1, pts_2);

    // get points 
    vector<pair<Vector3d, Vector3d>> corrs_3d, corrs_2d; 

    find3D2DMatch(pts_1, pts_2, depth_img1, corrs_3d, corrs_2d); 

    // cout <<"gt rotation: "<<endl<<rotation<<endl; 

    Matrix3d Rij = Matrix3d::Identity(); 
    Vector3d tij = Vector3d::Zero(); 

    // 3d-2d 
    MotionEstimator me;
    SolveTranslate st; 
    Matrix3d Rij_e, Rji_e, dR; //, Rij_e_t;
    Vector3d tij_e, tji_e, dt; // tij_e_t; 

    OptSolver opt_solver; 

    {
        vector<pair<Vector3d, Vector3d>> in_3d = getN(corrs_3d, cnt_3d); 
        vector<pair<Vector3d, Vector3d>> in_2d = in_3d;
        vector<pair<Vector3d, Vector3d>> in_2d_tmp;
        if(cnt_2d > cnt_3d){
            in_2d_tmp = getN(corrs_2d, cnt_2d-cnt_3d); 
            in_2d.insert(in_2d.end(), in_2d_tmp.begin(), in_2d_tmp.end()); 
        }

        // 3d-2d result 
        me.solvePNP_3D_2D(in_3d, Rij_e, tij_e); 
        print_err("opencv epnp: ", Rij_e, tij_e, Rij, tij); 

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

        cout<<"Hybrid PnP dt: "<<dt.norm()<<" dR: "<<computeAngle(dR)<<endl; 
        if(use_optimization){
            Matrix3d RR2 = Rij_e; 
            Vector3d tt2 = tij_e; 
            opt_solver.solveCeresHybrid(in_3d, in_2d_tmp, RR2, tt2); 
            dR = Rij.transpose()*RR2; 
            dt = tt2 - tij; 
            cout<<"Opti_2 dR: "<<dt.norm()<<" dR: "<<computeAngle(dR)<<endl; 
        }

        showMatchedFeatures(img_1, img_2, in_2d); 
        // v_dR.push_back(computeAngle(dR));  
        // v_dt.push_back(dt.norm()); 

    }
    return true; 
}

void showMatchedFeatures(cv::Mat img_1, cv::Mat img_2, vector<pair<Vector3d, Vector3d>>& matches)
{
    cv::Mat show_1, show_2; 
    cv::cvtColor(img_1, show_1, CV_GRAY2BGR); 
    cv::cvtColor(img_2, show_2, CV_GRAY2BGR); 
    std::vector< cv::KeyPoint > keypoints1;
    std::vector< cv::KeyPoint > keypoints2;
    vector<cv::DMatch> dmatches; 
    for(int i=0; i<matches.size(); i++){
        cv::KeyPoint pt1, pt2; 
        pt1.pt.x = matches[i].first(0) * FX + CX; 
        pt1.pt.y = matches[i].first(1) * FY + CY; 
        pt2.pt.x = matches[i].second(0) * FX + CX; 
        pt2.pt.y = matches[i].second(1) * FY + CY;
        keypoints1.push_back(pt1); 
        keypoints2.push_back(pt2); 
        cv::DMatch mm; 
        mm.queryIdx = mm.trainIdx = i; 
        dmatches.push_back(mm); 
    }

    cv::Mat show_img; 
    cv::drawMatches(show_1, keypoints1, show_2, keypoints2, dmatches, show_img); 

    cv::imshow("matches: ", show_img);
    cv::waitKey(0); 
}

void print_err(string pre, Matrix3d& Re, Vector3d& te, Matrix3d& Rg, Vector3d& tg)
{
	Vector3d dt = tg - te; 
	Matrix3d dR = Rg.transpose()*Re;

	cout<<pre<<" trans err: "<<dt.norm()<<" rotation err: "<<computeAngle(dR)<<endl; 
}