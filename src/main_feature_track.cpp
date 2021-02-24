/*
	Feb. 23, 2021, He Zhang, fuyinzh@gmail.com 

	Track features for several input images 

*/

#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <fstream>
#include <string>
#include <stdio.h>
#include <stdlib.h>

#include <opencv2/opencv.hpp>

using namespace Eigen; 
using namespace std; 

string folder("/home/davidz/work/data/iphone_egomotion"); 

// camera intrinsic params 
int COL = 1920; 
int ROW = 1440; 
double CX = 965.24;
double CY = 653.62; 
double FX = 1460.2;
double FY = 1460.2; 

// feature extraction related params 
double FEAT_MIN_DIST = 20; 
size_t MAX_NUM_FEAT = 2000; 
double FEAT_THRESHOLD = 0.01; // quality level 

double FUNDAMENTAL_MATRIX_THRESHOLD = 0.5; // used in rejectwithF 

void reduceVector(vector<cv::Point2f> &v, vector<uchar> status); 
bool inBorder(const cv::Point2f &pt);
vector<uchar> rejectWithF( vector<cv::Point2f>& pts_1, vector<cv::Point2f>& pts_2);
void showMatchedFeatures(cv::Mat img_1, cv::Mat img_2, vector<cv::Point2f>& pts_1, vector<cv::Point2f>& pts_2);

vector<vector<cv::Point2f>> feature_tracker(vector<string>& img_fnames, bool show_results = false); 
vector<double> get_depth(vector<cv::Point2f>& v_features, string dpt_file); 
bool save_to_file(string filename, const vector<vector<cv::Point2f>>& v_feats, const vector<double>& v_depths); 

int main(int argc, char* argv[])
{

    // folder = string("/media/davidz/Samsung_T5/dataset/gt_table/iphoneEgoMotion/rpy/x"); // y z
    folder = string("/media/davidz/Samsung_T5/dataset/gt_table/iphoneEgoMotion"); // y z

    vector<string> move_folders{"rpy/x/rot_roll", "rpy/y/rot_pitch", "rpy/z/rot_yaw", "trans_y/trans_y", "trans_z/trans_z"};
    vector<string> index{"1", "2", "3", "4", "5", "6", "7"};


    int image_num = 200; 
    int CNT = 100; 

    for(int k=0; k<move_folders.size(); k++){
        string move_folder = move_folders[k]; 
        for(int cnt=0; cnt < CNT; cnt++){
            vector<string> img_fnames; 
            for(int i=0; i<index.size(); i++){
                string img_file = folder + "/" + move_folder + "_" + index[i]+ "/color/" + to_string(image_num + cnt) + ".png"; 
                img_fnames.push_back(img_file); 
            }
            vector<vector<cv::Point2f>> features = feature_tracker(img_fnames, false); // true  
            string dpt_file = folder + "/" + move_folder + "_" + index[0]+ "/depth/" + to_string(image_num + cnt) + ".png"; 
            vector<double> depths = get_depth(features[0], dpt_file); 

            // save this 
            size_t found = move_folder.find_last_of("/"); 
            string file_out = "feats_tracker_results/" + move_folder.substr(found+1) + "_" + to_string(cnt+1); 
            if(!save_to_file(file_out, features, depths)){
                cerr<<" main_feature_track.cpp: failed to save file, stop!"; 
                return -1 ;
            }else{
                cout<<"succeed to save to "<<file_out<<endl; 
            }
        }
    }

    return 0; 
}

bool save_to_file(string filename, const vector<vector<cv::Point2f>>& v_feats, const vector<double>& v_depths)
{
    ofstream ouf(filename.c_str()); 
    if(!ouf.is_open()){
        cerr<<" main_feature_track.cpp: failed  to open file: "<<filename<<endl; 
        return false; 
    }
    // save file format, e.g. track features in 7 images, 1 2 3 4 5 6 7
        // depth (at 1st image), x1, y1, x2, y2, x3, y3, ..., x7, y7, normalized image planes 
    for(int i=0; i<v_depths.size(); i++){
        ouf<<v_depths[i]<<" "; 
        for(int j=0; j<v_feats.size(); j++){
            double xi = (v_feats[j][i].x - CX)/FX; 
            double yi = (v_feats[j][i].y - CY)/FY; 
            ouf<<xi<<" "<<yi<<" ";
        }
        ouf<<endl; 
    }
    ouf.flush();
    ouf.close(); 
    return true; 
}

vector<double> get_depth(vector<cv::Point2f>& v_features, string dpt_file)
{
    cv::Mat dpt_1 = cv::imread(dpt_file.c_str(), -1); 
    vector<double> v_depth(v_features.size(), -1.);
    double ui, vi, uj, vj; 
    double dpt_scale = 7.5; 
    double xi, yi, xj, yj;
    int cnt_3d = 0; 
    int cnt_2d = 0;  
    for(int i=0; i < v_features.size(); i++){

        ui = v_features[i].x; 
        vi = v_features[i].y; 

        int ui_dpt = round(ui/dpt_scale); 
        int vi_dpt = round(vi/dpt_scale); 

        double d = dpt_1.at<unsigned short>(vi_dpt, ui_dpt) * 0.001; 
        if(d >= 0.2 && d <= 7.5){
            // valid 3d point 
            v_depth[i] = d; 
            ++cnt_3d;
        }else
            ++cnt_2d; 
    }
    cout<<"main_feature_track.cpp: features: "<<v_features.size()<<" cnt_3d: "<<cnt_3d<<" cnt_2d: "<<cnt_2d<<endl; 
    return v_depth; 
}


vector<vector<cv::Point2f>> feature_tracker(vector<string>& img_fnames, bool show_track_result )
{
	std::vector<std::vector<cv::Point2f>> ret;
	vector<cv::Mat> img_v; // for display purpose 
	if(img_fnames.size() < 2){
		cerr<<" main_feature_track.cpp: input img_fnames.size(): "<<img_fnames.size(); 
		return ret; 
	}

	// extract the first feature pairs 
	vector<cv::Point2f> pts_1;
	cv::Mat img_1 = cv::imread(img_fnames[0].c_str(), -1); 
	cv::Mat gray_1; 
	cv::cvtColor(img_1, gray_1, CV_BGR2GRAY); 
	cv::goodFeaturesToTrack(gray_1, pts_1, MAX_NUM_FEAT, FEAT_THRESHOLD, FEAT_MIN_DIST); 

	vector<uchar> status;
    vector<float> err;
    vector<cv::Point2f> pts_2;
    cv::Mat img_2 = cv::imread(img_fnames[1].c_str(), -1); 
	cv::Mat gray_2; 
	cv::cvtColor(img_2, gray_2, CV_BGR2GRAY); 
    cv::calcOpticalFlowPyrLK(gray_1, gray_2, pts_1, pts_2, status, err, cv::Size(21, 21), 3);

    for (int i = 0; i < int(pts_2.size()); i++)
        if (status[i] && !inBorder(pts_2[i]))
                status[i] = 0;

    reduceVector(pts_1, status); 
    reduceVector(pts_2, status); 
    status = rejectWithF(pts_1, pts_2);
    reduceVector(pts_1, status); 
    reduceVector(pts_2, status); 

    // add to return tracked features 
    img_v.push_back(gray_1); 
    img_v.push_back(gray_2); 
    ret.push_back(pts_1); 
    ret.push_back(pts_2); 

    cv::Mat gray_cur = gray_2; 
    for(int j=2; j<img_fnames.size(); j++){
    	vector<cv::Point2f>& pts_cur = ret[j-1];
    	vector<cv::Point2f> pts_forw; 
    	vector<uchar> status;
    	cv::Mat img_forw = cv::imread(img_fnames[j].c_str(), -1);
    	cv::Mat gray_forw; 
    	cv::cvtColor(img_forw, gray_forw, CV_BGR2GRAY);
    	cv::calcOpticalFlowPyrLK(gray_cur, gray_forw, pts_cur, pts_forw, status, err, cv::Size(21, 21), 3);
    	for (int i = 0; i < int(pts_forw.size()); i++)
        	if (status[i] && !inBorder(pts_forw[i]))
                status[i] = 0;
        for(int k=0; k<j; k++)
        	reduceVector(ret[k], status);
        reduceVector(pts_forw, status);
        status = rejectWithF(pts_cur, pts_forw);
        for(int k=0; k<j; k++)
        	reduceVector(ret[k], status);
        reduceVector(pts_forw, status);
        ret.push_back(pts_forw); 
        gray_cur = gray_forw; 
        img_v.push_back(gray_forw); 
        cout <<" after tracking on image "<<j<<" feature matches number is: "<<ret[0].size()<<endl;
    }

    if(show_track_result){
        // see the result 
        for(int i=1; i<ret.size(); i++){
        	showMatchedFeatures(img_v[0], img_v[i], ret[0], ret[i]); 
        }
    }

    return ret; 

}

void showMatchedFeatures(cv::Mat img_1, cv::Mat img_2, vector<cv::Point2f>& pts_1, vector<cv::Point2f>& pts_2) 
{
    cv::Mat show_1, show_2; 
    cv::cvtColor(img_1, show_1, CV_GRAY2BGR); 
    cv::cvtColor(img_2, show_2, CV_GRAY2BGR);
    std::vector< cv::KeyPoint > keypoints1;
    std::vector< cv::KeyPoint > keypoints2;
    vector<cv::DMatch> dmatches; 
    for(int i=0; i<pts_1.size(); i++){
        cv::KeyPoint pt1, pt2; 
        pt1.pt = pts_1[i]; 
        pt2.pt = pts_2[i]; 
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


vector<uchar> rejectWithF( vector<cv::Point2f>& pts_1, vector<cv::Point2f>& pts_2)
{
	vector<uchar> status;
    if (pts_1.size() >= 8){
        double F_THRESHOLD = 0.5;
        cv::findFundamentalMat(pts_1, pts_2, cv::FM_RANSAC, F_THRESHOLD, 0.99, status);
        // reduceVector(pts_1, status);
        // reduceVector(pts_2, status);
    }else{
    	cerr<<" input feature matches are less than 8!"<<endl; 
    }

    return status; 
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

