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

double function(int k, double x){
}

double calculateMiddleRienmann(int K, vector<double> points, double(*fn)(int, double), int qn){
	double value = 0;
	for(int i = 1; i < points.size(); i++){
		value += (fn(K, (points(i-1) + points(i))/2) * (points(i) - points(i-1)));
	}
	return value;
}

double calculateTrapazoidalRienmann(int K, vector<double> points, double(*fn)(int, double), int qn){
	double value = 0;
	for(int i = 0; i < points.size(); i++){
		value += 0.5 * (fn(K, points(i-1)) + fn(K, points(i))) * (points(i) - points(i-1));
	}
	return value;
}

double calculateSimonRule(int K, vector<double> points, double(*fn)(int, double) , int qn){
	return ((2.0/3.0) * calculateMiddleRienmann(K,points,fn,qn) + (1.0/3.0) * calculateTrapazoidalRienmann(K,points,fn,qn));
}

int main(int argc, char *argv[]){

	return 0;
}
