/*

	Jan. 7, 2020, He Zhang, hzhang8@vcu.edu  

	some utility functions 

*/

#include "utility.h"

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