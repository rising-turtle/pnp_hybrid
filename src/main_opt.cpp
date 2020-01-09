/*
    Jan. 9th, 2020, He Zhang, hzhang8@vcu.edu 

    simulate for the optimization method with projected factors 

*/

#include "sim_corres.h"
#include "solve_5pts.h"
#include "solve_translate.h"
#include "utility.h"
#include <Eigen/Core>
#include <cmath>
#include <sstream>
#include <fstream>

using namespace Eigen; 
using namespace std; 

void 

int main(int argc, char* argv[])
{
    vector<int> v3d{4, 7, 10, 15, 20, 25, 30, 35, 40}; 
    run_monte_carlo(v3d, 40, 50); 

    return 0; 
}


