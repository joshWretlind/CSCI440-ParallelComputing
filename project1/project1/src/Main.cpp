#include<iostream>

using namespace std;

/**
 * Author: Josh Wretlind
 * Class: MATH/CSCI 440, Parallel Scientific Computing
 * Date: 09/13/2013
 * Assignment: 1) Serial Programming
 * Purpose: 
 **/

/**
 * createHilbertMatrix
 * Purpose: This is so that we can create a 2D array for a Hilbert matrix of N size
 * Arguments: int n: The size of the array that you would like to return
 * Return value: the 2D hilbert array of doubles
 **/
double** createHilbertMatrix(int n){
	//create a dynamic 2D array of N size
	double** hilbert = new double*[n];
	for(int i = 0; i < n; i++){
		hilbert[i] = new double[i];
	}

	//loop through the entire array, create the values that we need
	//equation for the values in a hibert matrix: A[i][j] = 1/(i + j - 1) i=1,2,...N, j=1,2,..N
	//Due to the first element in a array being indexed at 0, we must add 1 to i and j in the equation to get to their correct values
	//this forces the equation to be A[i][j] = 1/(i + j + 1) i = 0...N-1, j = 0...N-1
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++) {
			hilbert[i][j] = 1.0/(i + j + 1);
		}
	}

	return hilbert;
}

/**
 * This is main, the main entry point for an applicaiton.
 *
 **/
void main(){
	
}
