/*
	Jan. 28, 2021, He Zhang, fuyinzh@gmail.com
	
	wrapper to compute relative rotation given two sets of points  

*/


#include "rotation_only.h"

using namespace std;
using namespace Eigen;
using namespace opengv;

RotationOnly::RotationOnly(){}
RotationOnly::~RotationOnly(){}


namespace{

	void copy_to(vector<cv::Point2f>& fv, bearingVectors_t& tb){
		tb.resize(fv.size()); 
		for(int i=0; i<fv.size(); i++){
			Eigen::Vector3d v(fv[i].x, fv[i].y, 1.); 
			tb[i] = v; 
		}
	}

	void copy_to(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& matches, bearingVectors_t& ll, bearingVectors_t& rr){
		ll.resize(matches.size()); 
		rr.resize(matches.size()); 
		for(int i=0; i<matches.size(); i++){
			ll[i] = matches[i].first; 
			rr[i] = matches[i].second;
		}
	}

}

void RotationOnly::computeRotation(vector<cv::Point2f>& ll, std::vector<cv::Point2f>& rr, Eigen::Matrix3d& Rlr)
{
  	//derive correspondences based on random point-cloud
  	bearingVectors_t bearingVectors1;
  	bearingVectors_t bearingVectors2;

  	copy_to(ll, bearingVectors1); 
  	copy_to(rr, bearingVectors2); 
  	
  	//create a central relative adapter
  	relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2);

  	Rlr = relative_pose::rotationOnly(adapter);
  	return ; 
}



void RotationOnly::computeRotation(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& matches, Eigen::Matrix3d& Rlr)
{
  	//derive correspondences based on random point-cloud
  	bearingVectors_t bearingVectors1;
  	bearingVectors_t bearingVectors2;

  	copy_to(matches, bearingVectors1, bearingVectors2);
  	
  	//create a central relative adapter
  	relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2);

  	Rlr = relative_pose::rotationOnly(adapter);
  	return ; 
}

void RotationOnly::computeRotationSac(vector<cv::Point2f>& ll, std::vector<cv::Point2f>& rr, Eigen::Matrix3d& Rlr)
{
	//derive correspondences based on random point-cloud
  	bearingVectors_t bearingVectors1;
  	bearingVectors_t bearingVectors2;

  	copy_to(ll, bearingVectors1); 
  	copy_to(rr, bearingVectors2); 

  	  //create a central relative adapter
  	relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2);

  	//Create a RotationOnlySacProblem and Ransac
	sac::Ransac<
	sac_problems::relative_pose::RotationOnlySacProblem> ransac;
	std::shared_ptr<
	sac_problems::relative_pose::RotationOnlySacProblem> relposeproblem_ptr(
	new sac_problems::relative_pose::RotationOnlySacProblem(adapter));
	ransac.sac_model_ = relposeproblem_ptr;
	ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
	ransac.max_iterations_ = 50;
	ransac.computeModel(0);
	// std::cout << ransac.model_coefficients_ << std::endl << std::endl;
	Rlr = ransac.model_coefficients_; 
	return ; 
}
void RotationOnly::computeRotationSac(std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>>& matches, Eigen::Matrix3d& Rlr)
{
	//derive correspondences based on random point-cloud
  	bearingVectors_t bearingVectors1;
  	bearingVectors_t bearingVectors2;

	copy_to(matches, bearingVectors1, bearingVectors2);

  	  //create a central relative adapter
  	relative_pose::CentralRelativeAdapter adapter(
      bearingVectors1,
      bearingVectors2);

  	//Create a RotationOnlySacProblem and Ransac
	sac::Ransac<
	sac_problems::relative_pose::RotationOnlySacProblem> ransac;
	std::shared_ptr<
	sac_problems::relative_pose::RotationOnlySacProblem> relposeproblem_ptr(
	new sac_problems::relative_pose::RotationOnlySacProblem(adapter));
	ransac.sac_model_ = relposeproblem_ptr;
	ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
	ransac.max_iterations_ = 50;
	ransac.computeModel(0);
	// std::cout << ransac.model_coefficients_ << std::endl << std::endl;
	Rlr = ransac.model_coefficients_; 
	return ; 
}