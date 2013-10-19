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
int totalSize;
int kPerWorker;
int lowerKBound;
int upperBound;
long double regionPerK;


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
long double function(int k, long double x){
	return cos(100.0 * x - k * sin(x));
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
long double calculateMiddleRienmann(int K, vector<long double> points, long double(*fn)(int, long double), int qn){
	double value = 0;
	for(int i = 1; i < points.size(); i++){
		value += (fn(K, (points[i-1] + points[i])*0.5) * (points[i] - points[i-1]));
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
long double calculateTrapazoidalRienmann(int K, vector<long double> points, long double(*fn)(int, long double), int qn){
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
long double calculateSimonRule(int K, vector<long double> points, long double(*fn)(int, long double) , int qn){
	return ((2.0/3.0) * calculateMiddleRienmann(K,points,fn,qn) + (1.0/3.0) * calculateTrapazoidalRienmann(K,points,fn,qn));
}

vector<long double> providePoints(long double begin, long double end, int regions){
	vector<long double> points;
	points.resize(regions + 1);
	for(int i = 0; i <= regions; i++){
		points[i] = begin + i*(end - begin)/((long double)regions);
	}
	
	return points;
}

long double calculateMiddleSum(int k, int rank){
	double sum = 0;
	for(int j = lowerKBound; j < upperBound; j++){
		vector<long double> points = providePoints(j*regionPerK,(j+1)*regionPerK,100);
		sum += calculateMiddleRienmann(k,points,function,points.size());
	}
	return sum;
}

long double calculateSimonSum(int k, int rank){
	long double sum = 0;
	for(int j = lowerKBound; j < upperBound; j++){
		vector<long double> points = providePoints(j*regionPerK,(j+1)*regionPerK,100);
		sum += calculateSimonRule(k,points,function,points.size());
	}
	return sum;
}

long double calculateTrapazoidSum(int k, int rank){
	long double sum = 0;
	for(int j = lowerKBound; j < upperBound; j++){
		vector<long double> points = providePoints(j*regionPerK,(j+1)*regionPerK,100);
		sum += calculateTrapazoidalRienmann(k,points,function,points.size());
	}
	return sum;
}

long double* calculateIntegral(int k, int rank){
	kPerWorker = ceil(((long double)k)/((long double)totalSize));
	lowerKBound = rank*(kPerWorker);
	upperBound = (rank+1)*(kPerWorker);
	if((rank+1) == totalSize){
		upperBound = k;
	}
	if(upperBound > k){
		upperBound = k;
	}
	if(lowerKBound > k){
		lowerKBound = k;
	}
	if(lowerKBound == upperBound){
		static long double arr[3] = {0,0,0};
		return arr;
	}
	regionPerK = M_PI/((long double)k);
	cout << "My Rank: " << rank << " LowerBound: " << lowerKBound << " Upper: " << upperBound << " perK: " << regionPerK << endl;
	static long double results[3] = {0,0,0};
	results[0] = calculateMiddleSum(k,rank);
	results[1] = calculateTrapazoidSum(k,rank);
	results[2] = calculateSimonSum(k,rank);
	return results;
}

int main(int argc, char *argv[]){
	MPI::Init(argc, argv);
	
	totalSize = MPI::COMM_WORLD.Get_size();
	int myRank = MPI::COMM_WORLD.Get_rank();
	long double *k = calculateIntegral(100,myRank);
	long double sum[3] = {0,0,0};
	MPI::COMM_WORLD.Reduce(&k,&sum,3,MPI::LONG_DOUBLE,MPI_SUM,master);
	
	if(myRank == 0){
		cout << "Reduced count: " << sum[0] << " " << sum[1] << " " << sum[2] << endl;
	}
	MPI::Finalize();
	return 0;
}
