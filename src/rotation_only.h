/*
	Jan. 28, 2021, He Zhang, fuyinzh@gmail.com
	
	wrapper to compute relative rotation given two sets of points  

*/

#pragma once 


#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <limits.h>
#include <vector>
#include <map>
#include <Eigen/Eigen>
#include <opencv2/opencv.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <opengv/sac/Ransac.hpp>
#include <opengv/sac/Lmeds.hpp>
#include <opengv/sac_problems/relative_pose/RotationOnlySacProblem.hpp>
#include <sstream>
#include <fstream>


class RotationOnly{

public:
	RotationOnly(); 
	~RotationOnly(); 


	void computeRotation( std::vector<cv::Point2f>& ll, std::vector<cv::Point2f>& rr, Eigen::Matrix3d& Rlr);
	void computeRotation( std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d> >& matches, Eigen::Matrix3d& Rlr);

	void computeRotationSac( std::vector<cv::Point2f>& ll, std::vector<cv::Point2f>& rr, Eigen::Matrix3d& Rlr);
	void computeRotationSac( std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d> >& matches, Eigen::Matrix3d& Rlr);


};