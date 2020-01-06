/*

    Jan. 6th, 2020, He Zhang, hzhang8@vcu.edu 
    
    simulate point correspondences given [R, t]

*/

#include "sim_corres.h"
#include <random>

SimCorr::SimCorr(){
    init(); 
}
SimCorr::~SimCorr(){}


void SimCorr::init(int pix_step)
{
    // assume fx = 460, cx = 320, cy = 240, 
    float cx = 320; 
    float cy = 240; 
    float f = 460; 
    static std::default_random_engine gen;
    
    // init feature points randomly locates in range z-axis [2-7] meters 
    // Initializing of uniform_real_distribution class 
    uniform_real_distribution<double> uniform_dis(2, 7); 
    
    mg_feats.clear(); 
    int id = 0; 
    for(int r=pix_step; r< 480; r+=pix_step)
    for(int c=pix_step; c< 640; c+=pix_step)
    {
	Vector3d p( (c-cx)/f,  (r-cy)/f, 1.); 
	p.z() = uniform_dis(gen); 
	p.x() *= p.z(); 
	p.y() *= p.z(); 
	mg_feats.emplace(make_pair(++id, p));
    }
    
    cout<<"sim_corres.cpp: init "<<id<<" features "<<endl; 

    return ; 
}

vector<pair<Vector3d, Vector3d>> SimCorr::find_corrs( Matrix3d& Rji, Vector3d& tji)
{
    vector<pair<Vector3d, Vector3d>> ret; 
    map<int, Vector3d>::iterator it = mg_feats.begin();
    while(it != mg_feats.end()){
	
	Vector3d pi = it->second; 
	Vector3d pj = Rji*pi + tji; 
	
	if(pj.z() >= 0.2){
	    
	    ret.push_back(make_pair(pi, pj)); 
	} 
	++it; 
    }
    
    return ret; 
}


// add noise to the corresponse features 
vector<pair<Vector3d, Vector3d>> SimCorr::add_noise( vector<pair<Vector3d, Vector3d>>& corres, double pix_std, poly d_std)
{
    vector<pair<Vector3d, Vector3d>> ret; 
    static std::default_random_engine gen;
    std::normal_distribution<double> pix_gauss(0.0,pix_std);

    for(int i=0; i<corres.size(); i++){

	Vector3d pi = corres[i].first; 
	Vector3d pj = corres[i].second; 
	
	pi.x() += pix_gauss(gen); 
	pi.y() += pix_gauss(gen); 
	
	pj.x() += pix_gauss(gen); 
	pj.y() += pix_gauss(gen); 
	
	double dpt_i = pi.z(); 
	double dpt_j = pj.z(); 
	

	std::normal_distribution<double> dpt_i_gauss(0.0, d_std.y(dpt_i)); 
	std::normal_distribution<double> dpt_j_gauss(0.0, d_std.y(dpt_j)); 
	
	pi.z() += dpt_i_gauss(gen); 
	pj.z() += dpt_j_gauss(gen); 

	ret.push_back(make_pair(pi, pj)); 
    }

    return ret; 

}




