/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu 

	solve t in [R,t] given R 
	
*/

#include "solve_translate.h"
#include "translate_factor.h"
#include "utility.h"

namespace{

	// based on the 
	void compute_cov(Matrix3d& cov, Vector3d& pt){

		static poly std_poly; 

		double sig_z = std_poly.y(pt.z()); 

		Vector3d np(pt.x(), pt.y(), 1); 
		cov = SQ(sig_z)*np*np.transpose();

		double sig_u = (1./460); // 1 pixel deviation
		double sig_v = (1./460); 

		cov(0,0) += SQ(pt.z())*SQ(sig_u); 
		cov(1,1) += SQ(pt.z())*SQ(sig_v); 
		return ; 
	}

}

SolveTranslate::SolveTranslate(){}

bool SolveTranslate::solveTCeres(const vector<pair<Vector3d, Vector3d>> &corres, const Matrix3d &Rij, Vector3d &tij)
{
    Matrix3d Rji = Rij.transpose(); 
	Eigen::Vector3d tji = -Rji*tij; 
	para_T[0][0] = tji(0); 
	para_T[0][1] = tji(1); 
	para_T[0][2] = tji(2); 

  	ceres::Problem problem;
    ceres::LossFunction *loss_function;
    // loss_function = new ceres::CauchyLoss(1.0);    
    loss_function = new ceres::HuberLoss(1.0);

    problem.AddParameterBlock(para_T[0], 3); 

    Matrix3d cov_pti = Matrix3d::Identity(); 

    for(int i=0; i<corres.size(); i++){

    	Vector3d pti = corres[i].first; 
    	Vector3d ptj = corres[i].second; 

    	compute_cov(cov_pti, pti); 

    	pti.x() *= pti.z(); 
    	pti.y() *= pti.z();
    	Vector2d ptj_n(ptj(0), ptj(1)); 

    	TranslateFactor *f = new TranslateFactor(Rji, pti, ptj_n, cov_pti); 

    	ceres::ResidualBlockId fid = problem.AddResidualBlock(f, loss_function, para_T[0]); 

    	// if(i< 10){
    	// 	 vector<double*>* para = new vector<double*>; 
    	// 	 problem.GetParameterBlocksForResidualBlock(fid, para); 
     //         vector<double> res(2); 
     //         f->Evaluate(&para[0][0], &res[0], 0); 
     //         cout<<"solve_translate.cpp: residual: "<<res[0]<<" "<<res[1]<<endl;
    	// }
    }

    ceres::Solver::Options options;
    options.linear_solver_type = ceres::DENSE_SCHUR;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = 20;
    // options.minimizer_progress_to_stdout = true;
    // options.max_solver_time_in_seconds = SOLVER_TIME; 
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    // cout << summary.BriefReport() << endl;

    // Vector3d tji(para_T[0][0], para_T[0][1], para_T[0][2]); 
    tji = Vector3d(para_T[0][0], para_T[0][1], para_T[0][2]);
    tij = -Rij*tji;

	return true ; 

}