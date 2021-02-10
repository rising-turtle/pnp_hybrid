/*
	Jan. 9, 2020, He Zhang, hzhang8@vcu.edu 

	functions related to sfm 

*/


#include "opt_solver.h"
#include "utility.h"

OptSolver::OptSolver(){}


	// combine 2D and 3D 
bool OptSolver::solveCeresHybrid(const vector<pair<Vector3d, Vector3d>> &corres_3d, 
		const vector<pair<Vector3d, Vector3d>> &corres_2d, Matrix3d &Rij, Vector3d &tij)
{
	cout<<"before opt: tij: "<<tij.transpose()<<endl; 
	cout<<"sum 3D err: "<<sum_error(corres_3d, Rij, tij)<<" 2D err: "<<sum_error_2d(corres_2d, Rij, tij)<<endl; 

	Matrix3d Rji = Rij.transpose(); 
	Eigen::Vector3d tji = -Rji*tij; 
	Quaterniond qji(Rji); 

	m_tji[0][0] = tji.x(); 
	m_tji[0][1] = tji.y(); 
	m_tji[0][2] = tji.z(); 

	m_qji[0][0] = qji.w(); 
	m_qji[0][1] = qji.x(); 
	m_qji[0][2] = qji.y(); 
	m_qji[0][3] = qji.z(); 
 	
 
	//full BA
	ceres::Problem problem;
	ceres::LocalParameterization* local_parameterization = new ceres::QuaternionParameterization();
	problem.AddParameterBlock(m_qji[0], 4, local_parameterization);
	problem.AddParameterBlock(m_tji[0], 3);

	int N = corres_3d.size(); 
	for(int i=0; i<N; i++){
		Vector3d pti = corres_3d[i].first; 
		pti.x() *= pti.z(); 
    	pti.y() *= pti.z();
    	para_pt[i][0] = pti.x(); 
    	para_pt[i][1] = pti.y(); 
    	para_pt[i][2] = pti.z(); 
    	problem.AddParameterBlock(para_pt[i], 3);
    	problem.SetParameterBlockConstant(para_pt[i]);
	}


	for(int i=0; i<corres_3d.size(); i++){


		ceres::CostFunction* cost_function = ReprojectionError3D::Create(
											corres_3d[i].second.x(),
											corres_3d[i].second.y());
		problem.AddResidualBlock(cost_function, NULL, m_qji[0], m_tji[0], para_pt[i]); 

	}

	for(int i=0; i<corres_2d.size(); i++){
		ceres::CostFunction* cost_function = SampsonEssential::Create(
											corres_2d[i].first.x(),
											corres_2d[i].first.y(),
											corres_2d[i].second.x(),
											corres_2d[i].second.y());
		problem.AddResidualBlock(cost_function, NULL, m_qji[0], m_tji[0]); 
	}

	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.max_solver_time_in_seconds = 0.2;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	// cout << summary.BriefReport() << endl;

	qji.w() = m_qji[0][0]; 
	qji.x() = m_qji[0][1]; 
	qji.y() = m_qji[0][2]; 
	qji.z() = m_qji[0][3]; 

	Quaterniond qij = qji.inverse(); 
	tji.x() = m_tji[0][0]; 
	tji.y() = m_tji[0][1]; 
	tji.z() = m_tji[0][2]; 
	tij = qij * tji;
	tij*=-1.;

	Rij = qij.toRotationMatrix();

	cout<<"after opt: tij: "<<tij.transpose()<<endl; 
	cout<<"sum 3D err: "<<sum_error(corres_3d, Rij, tij)<<" 2D err: "<<sum_error_2d(corres_2d, Rij, tij)<<endl; 

	return true; 
}

bool OptSolver::solveCeres(const vector<pair<Vector3d, Vector3d>> &corres, Matrix3d &Rij, Vector3d &tij)
{

	cout<<"before opt: tij: "<<tij.transpose()<<endl; 
	cout<<"sum err: "<<sum_error(corres, Rij, tij)<<endl; 

	Matrix3d Rji = Rij.transpose(); 
	Eigen::Vector3d tji = -Rji*tij; 
	Quaterniond qji(Rji); 

	m_tji[0][0] = tji.x(); 
	m_tji[0][1] = tji.y(); 
	m_tji[0][2] = tji.z(); 

	m_qji[0][0] = qji.w(); 
	m_qji[0][1] = qji.x(); 
	m_qji[0][2] = qji.y(); 
	m_qji[0][3] = qji.z(); 
 	
 
	//full BA
	ceres::Problem problem;
	ceres::LocalParameterization* local_parameterization = new ceres::QuaternionParameterization();
	problem.AddParameterBlock(m_qji[0], 4, local_parameterization);
	problem.AddParameterBlock(m_tji[0], 3);

	int N = corres.size(); 
	for(int i=0; i<N; i++){
		Vector3d pti = corres[i].first; 
		pti.x() *= pti.z(); 
    	pti.y() *= pti.z();
    	para_pt[i][0] = pti.x(); 
    	para_pt[i][1] = pti.y(); 
    	para_pt[i][2] = pti.z(); 
    	problem.AddParameterBlock(para_pt[i], 3);
    	problem.SetParameterBlockConstant(para_pt[i]);
	}


	for(int i=0; i<corres.size(); i++){


		ceres::CostFunction* cost_function = ReprojectionError3D::Create(
											corres[i].second.x(),
											corres[i].second.y());
		problem.AddResidualBlock(cost_function, NULL, m_qji[0], m_tji[0], para_pt[i]); 

	}

	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;
	options.max_solver_time_in_seconds = 0.2;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	// cout << summary.BriefReport() << endl;

	qji.w() = m_qji[0][0]; 
	qji.x() = m_qji[0][1]; 
	qji.y() = m_qji[0][2]; 
	qji.z() = m_qji[0][3]; 

	Quaterniond qij = qji.inverse(); 
	tji.x() = m_tji[0][0]; 
	tji.y() = m_tji[0][1]; 
	tji.z() = m_tji[0][2]; 
	tij = qij * tji;
	tij*=-1.;

	Rij = qij.toRotationMatrix();

	cout<<"after opt: tij: "<<tij.transpose()<<endl; 
	cout<<"sum err: "<<sum_error(corres, Rij, tij)<<endl; 

	return true; 
}
