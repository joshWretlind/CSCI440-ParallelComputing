/*************************************
 * Project 3
 * Purpose: This is to answer the second part of the assignment
 * Author: Josh Wretlind
 * Date: 10/17/2013
 *************************************/
#include<iostream>
#include "mpi.h"
#include<string>
#include<vector>
#include<cmath>

using namespace std;
using namespace MPI;

/*************************************
 * function(int k, vector<double> points)
 * Purpose: This is to calculate a collection of values corrosponding to points
 *          action function it calculates: f(x) = cos(100x - k*sin(x))
 * Parameters: int k: The wave number for the function we want to calculate
 *             vector<double> points: the collection of points we want to evaluate the function on
 * Output: vector<double> of the values corrosponsing with points
 * Complexity: Time:  O(N)
 *             Space: O(N) 
 ***********************************/
vector<double> function(int k, vector<double> points){
	vector<double> results;
	results.resize(points.size());
	for(int i = 0; i < points.size(); i++){
		results[i] = cos(100 * points[i] - k * sin(points[i]));
	}
	return results;
}

int main(int argc, char *argv[]){

	return 0;
}
