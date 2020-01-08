/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu  

	some utility functions 

*/

#include "utility.h"
#include <cmath>
#include <algorithm>

using namespace std; 

vector<pair<Vector3d, Vector3d>> getN(const vector<pair<Vector3d, Vector3d>>& in, int N){
	vector<pair<Vector3d, Vector3d>> ret; 

	int tN = in.size(); 
	int step = (tN/N);
	if(step <= 0) step = 1; 

	for(int i=0; i<tN; i+=step){
		ret.push_back(make_pair(in[i].first, in[i].second)); 
	}

	return ret; 
}

double computeAngle(Matrix3d& R)
{
	double dtrace = (R(0,0) + R(1,1) + R(2,2)-1.)/2.; 
	if(dtrace< -1) dtrace = -1; 
	if(dtrace>1) dtrace = 1; 
	return acos(dtrace);
}

pair<double, double> getMeanStd(vector<double>& v)
{
	double sum = std::accumulate(std::begin(v), std::end(v), 0.0);
	double m =  sum / v.size();

	double accum = 0.0;
	std::for_each (std::begin(v), std::end(v), [&](const double d) {
	    accum += (d - m) * (d - m);
	});

	double stdev = sqrt(accum / (v.size()-1));
	return make_pair(m, stdev); 
}

vector<pair<Vector3d, Vector3d>> getInliers(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask)
{
    assert(mask.rows == corres.size()); 

    vector<pair<Vector3d, Vector3d>> ret;
    for(int i=0; i<corres.size(); i++){

        if(mask.at<unsigned char>(i,0) == 0) continue; 
        
        ret.push_back(make_pair(corres[i].first, corres[i].second)); 
    } 
    return ret; 
}

vector<pair<Vector3d, Vector3d>> getInliersIndex(const vector<pair<Vector3d, Vector3d>> &corres, cv::Mat& mask)
{
	vector<pair<Vector3d, Vector3d>> ret;
	int index; 
    for(int i=0; i<mask.rows; i++){

    	index = mask.at<unsigned char>(i, 0); 

    	if(index > corres.size()){
    		cerr<<"utility.cpp: index: "<<index<<" > corres.size(): "<<corres.size()<<endl; 
    		continue; 
    	}
        
        ret.push_back(make_pair(corres[index].first, corres[index].second)); 
    } 
    return ret; 
}     
