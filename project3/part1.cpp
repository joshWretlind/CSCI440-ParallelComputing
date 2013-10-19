/*************************************
 * Project 3
 * Purpose: This is to answer the first part of the assignment
 * Author: Josh Wretlind
 * Date: 10/17/2013
 *************************************/
#include<iostream>
#include "mpi.h"
#include<string>
#include<vector>

using namespace std;
using namespace MPI;

/*************************************
 * function(int k, double x)
 * Purpose: This is the dummy function declaration for this part of the project.
 *          Basically, what is supposed to happen is that you call the calculate methods
 *          with this name in the fn parameter. Ex: calculateSimonRule(100,[1,2,3,4], &function, 4)
 * Parameters: int k: This is the wave number
 *             int x: this is the x value that the function should take
 * Output: The value of the fucntion with a wave number of k at x
 * Complexity: Time:  Unknown at the moment, this function hasn't been filled in, try for O(1)
 *             Space: Unknwon, try for O(1)
 *************************************/
double function(int k, double x){
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

	return 0;
}
