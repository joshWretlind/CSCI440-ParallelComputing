/*************************************
 * Project 3
 * Purpose: This is to answer the third part of the assignment
 * Author: Josh Wretlind
 * Date: 10/17/2013
 *************************************/
#include<iostream>
#include "mpi.h"
#include<string>
#include<cmath> 
#include<vector>

using namespace std;

int master = 0;

/*************************************
 * function(int k, double points)
 * Purpose: This is to calculate a  values corrosponding to a given point
 *          function it calculates: f(x) = cos(100x - k*sin(x))
 * Parameters: int k: The wave number for the function we want to calculate
 *             double points: the  points we want to evaluate the function for
 * Output: double of the value corrosponsing with the point given
 * Complexity: Time:  O(N)
 *             Space: O(N) 
 *************************************/
double function(int k, double x){
	return cos(100 * x - k * sin(x));
}

/*************************************
 * calculateMiddleRienmann(int K, vector<double> points, function* fn, int qn)
 * Purpose: This method calculates the Rienmann middle sum for fn using the midpoints
 *          of points and a wave number of K.
 * Parameters: int K: This is the wave number for the function we're calculating the sum for
 *             vector<double> points: this is the collection of points that we should use for calculating the sum
 *             double(*fn)(int, double): This is the pointer to the method to calculate 
 *                                       the values for the function we want to integrate
 *             int qn: This is the size of the number of points we want to calculate
 * Output: The calculated Rienmann middle sum
 * Complexity: Time:  O(N) where N is the size of points
 *             Space: O(1)
 *************************************/
double calculateMiddleRienmann(int K, vector<double> points, double(*fn)(int, double), int qn){
	double value = 0;
	for(int i = 1; i < points.size(); i++){
		value += (fn(K, (points[i-1] + points[i])/2) * (points[i] - points[i-1]));
	}
	return value;
}

/*************************************
 * calculateTrapazoidalRienmann(int K, vector<double> points, function* fn, int qn)
 * Purpose: This method calculates the trapazoidal sum for fn using the midpoints
 *          of points and a wave number of K.
 * Parameters: int K: This is the wave number for the function we're calculating the sum for
 *             vector<double> points: this is the collection of points that we should use for calculating the sum
 *             double(*fn)(int, double): This is the pointer to the method to calculate 
 *                                       the values for the function we want to integrate
 *             int qn: This is the size of the number of points we want to calculate
 * Output: The calculated trapazoidal sum
 * Complexity: Time:  O(N) where N is the size of points
 *             Space: O(1)
 *************************************/
double calculateTrapazoidalRienmann(int K, vector<double> points, double(*fn)(int, double), int qn){
	double value = 0;
	for(int i = 0; i < points.size(); i++){
		value += 0.5 * (fn(K, points[i-1]) + fn(K, points[i])) * (points[i] - points[i-1]);
	}
	return value;
}

/*************************************
 * calculateSimonRule(int K, vector<double> points, function* fn, int qn)
 * Purpose: This method calculates the first order simon rule for fn using the midpoints
 *          of points and a wave number of K.
 * Parameters: int K: This is the wave number for the function we're calculating the sum for
 *             vector<double> points: this is the collection of points that we should use for calculating the sum
 *             double(*fn)(int, double): This is the pointer to the method to calculate 
 *                                       the values for the function we want to integrate
 *             int qn: This is the size of the number of points we want to calculate
 * Output: The calculated first order simon rule
 * Complexity: Time:  O(N) where N is the size of points
 *             Space: O(1)
 *************************************/
double calculateSimonRule(int K, vector<double> points, double(*fn)(int, double) , int qn){
	return ((2.0/3.0) * calculateMiddleRienmann(K,points,fn,qn) + (1.0/3.0) * calculateTrapazoidalRienmann(K,points,fn,qn));
}


int main(int argc, char *argv[]){
	MPI::Init(argc, argv);
	
	int totalSize = MPI::COMM_WORLD.Get_size();
	int myRank = MPI::COMM_WORLD.Get_rank();
	int k = 100;
	int sum = 0;
	MPI::COMM_WORLD.Reduce(k,sum,1,MPI::INTEGER,MPI_SUM);
		
	MPI::Finalize();
	return 0;
}
