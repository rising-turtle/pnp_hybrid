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
	// p.x() *= p.z(); 
	// p.y() *= p.z(); 
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
	pi.x() = pi.x()*pi.z(); 
	pi.y() = pi.y()*pi.z();  
	Vector3d pj = Rji*pi + tji; 
	
	if(pj.z() >= 0.2){
	    	
	    pj.x() = pj.x() / pj.z(); 
	    pj.y() = pj.y() / pj.z(); 

	    ret.push_back(make_pair(it->second, pj)); 
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

// add noise to the corresponse features 
vector<pair<Vector3d, Vector3d>> SimCorr::addNoise3D2D( vector<pair<Vector3d, Vector3d>>& corres, double pix_std, poly d_std)
{
    vector<pair<Vector3d, Vector3d>> ret; 
    static std::default_random_engine gen;
    std::normal_distribution<double> pix_gauss(0.0,pix_std);

    for(int i=0; i<corres.size(); i++){

		Vector3d pi = corres[i].first; 
		Vector3d pj = corres[i].second; 

		std::normal_distribution<double> dpt_i_gauss(0.0, d_std.y(pi.z())); 
		double z = pi.z(); //  + dpt_i_gauss(gen); 

		pi = pi/pi.norm(); 
		pj = pj/pj.norm(); 

		pi = addNoiseRay(pix_std, pi); 
		pj = addNoiseRay(pix_std, pj); 
		
		pi = (pi/pi.z()) * z; 
		
		ret.push_back(make_pair(pi, pj)); 
    }

    return ret; 

}


Eigen::Vector3d SimCorr::addNoiseRay( double noiseLevel, Eigen::Vector3d cleanPoint )
{
  //compute a vector in the normal plane (based on good conditioning)
  Eigen::Vector3d normalVector1;
  if(
      (fabs(cleanPoint[0]) > fabs(cleanPoint[1])) &&
      (fabs(cleanPoint[0]) > fabs(cleanPoint[2])) )
  {
    normalVector1[1] = 1.0;
    normalVector1[2] = 0.0;
    normalVector1[0] = -cleanPoint[1]/cleanPoint[0];
  }
  else
  {
    if(
        (fabs(cleanPoint[1]) > fabs(cleanPoint[0])) &&
        (fabs(cleanPoint[1]) > fabs(cleanPoint[2])) )
    {
      normalVector1[2] = 1.0;
      normalVector1[0] = 0.0;
      normalVector1[1] = -cleanPoint[2]/cleanPoint[1];
    }
    else
    {
      normalVector1[0] = 1.0;
      normalVector1[1] = 0.0;
      normalVector1[2] = -cleanPoint[0]/cleanPoint[2];
    }
  }

  normalVector1 = normalVector1 / normalVector1.norm();
  Eigen::Vector3d normalVector2 = cleanPoint.cross(normalVector1);
  double noiseX =
      noiseLevel * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;
  double noiseY =
      noiseLevel * (((double) rand())/ ((double) RAND_MAX)-0.5)*2.0 / 1.4142;

  // Eigen::Vector3d noisyPoint =
  //    800 * cleanPoint + noiseX *normalVector1 + noiseY * normalVector2;
  Eigen::Vector3d noisyPoint =
       cleanPoint + noiseX *normalVector1 + noiseY * normalVector2;
  noisyPoint = noisyPoint / noisyPoint.norm();
  return noisyPoint;

}


