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


double sum_error(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d& Rij, const Vector3d& tij)
{
	double ret = 0; 
	Matrix3d Rji = Rij.transpose(); 
	Vector3d tji = -Rji*tij; 

	// cout<<"Rij: "<<endl<<Rij<<endl<<tij.transpose()<<endl; 

	int N = corres.size(); 
	// Eigen::MatrixXd long_err(N*2, 1); 
	for(int i=0; i<corres.size(); i++){

		Vector3d pti = corres[i].first; 
		Vector3d ptj = corres[i].second; 

 		pti.x() *= pti.z(); 
    	pti.y() *= pti.z();

    	Vector3d pti_j = Rji*pti + tji; 
    	pti_j.x() /= pti_j.z(); 
    	pti_j.y() /= pti_j.z(); 
    	Vector2d err(pti_j.x() - ptj.x(), pti_j.y() - ptj.y()); 
    	// long_err(i*2) = err(0); 
    	// long_err(i*2+1) = err(1); 
    	// err = err * pti_j.z();
    	ret += err.squaredNorm(); 
    	// cout <<" in sum_error i: "<<i<<" err: "<<err.transpose()<<" err.norm: "<<err.norm()<<endl; 
	}
	return 0.5*ret; 
	// return 0.5*long_err.squaredNorm(); 
}

template <typename Derived>
static Eigen::Matrix<typename Derived::Scalar, 3, 3> skewSymmetric(const Eigen::MatrixBase<Derived> &q)
{
    Eigen::Matrix<typename Derived::Scalar, 3, 3> ans;
    ans << typename Derived::Scalar(0), -q(2), q(1),
        q(2), typename Derived::Scalar(0), -q(0),
        -q(1), q(0), typename Derived::Scalar(0);
    return ans;
}

double sum_error_2d(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d& Rij, const Vector3d& tij)
{
	double ret = 0; 
	Matrix3d Rji = Rij.transpose(); 
	Vector3d tji = -Rji*tij; 

	Matrix3d E = Rji * skewSymmetric(tij);
	//cout<<"Rji: "<<endl<<Rji<<endl;
	//cout<<"tij: "<<endl<<tij<<endl;

	// cout<<"E: "<<endl<<E<<endl; 

	int N = corres.size(); 
	for(int i=0; i<corres.size(); i++){

		Vector3d pti = corres[i].first; 
		Vector3d ptj = corres[i].second; 

		pti(2) = 1.; 
		ptj(2) = 1.; 

		// cout<<"pti: "<<pti.transpose()<<" ptj: "<<ptj.transpose()<<endl; 

		double epsilon = ptj.transpose() * E * pti; 
		Eigen::Vector3d ep_i = E * pti; 
    	Eigen::Vector3d ep_j = E.transpose() * ptj; 

    	//cout<<"ep_i: "<<ep_i.transpose()<<" ep_j: "<<ep_j.transpose()<<endl; 

    	Eigen::Matrix<double, 1, 4> J; 
    	J << ep_j(0), ep_j(1), ep_i(0), ep_i(1);
    	double JJ = ep_i(0)*ep_i(0) + ep_i(1)*ep_i(1) + ep_j(0)*ep_j(0) + ep_j(1)*ep_j(1); 
    	
    	// sampson approximation 
    	Eigen::Matrix<double, 4, 1> residual;
    	double inv_JJ = 1./JJ; 
    	residual = -J.transpose() * inv_JJ * epsilon; 
    	ret += residual.squaredNorm(); 
	}
	return 0.5*ret; 
}