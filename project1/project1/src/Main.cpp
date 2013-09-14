#include<iostream>
#include<cmath>

using namespace std;

/**
 * Author: Josh Wretlind
 * Class: MATH/CSCI 440, Parallel Scientific Computing
 * Date: 09/13/2013
 * Assignment: 1) Serial Programming
 * Purpose: This calculates a few matricies, as well as does matrix-matrix multiplicaiton on them.
 **/

/**
 * createHilbertMatrix
 * Purpose: This is so that we can create a 2D array for a Hilbert matrix of N size
 * Arguments: int n: The size of the array that you would like to return
 * Return value: the 2D hilbert array of doubles
 * Complexity: Time:  O(N^2)
 *             Space: O(N^2)
 **/
long double** createHilbertMatrix(int n){
	//create a dynamic 2D array of N size
	long double** hilbert = new long double*[n];
	for(int i = 0; i < n; i++){
		hilbert[i] = new long double[i];
	}

	//loop through the entire array, create the values that we need
	//equation for the values in a hibert matrix: A[i][j] = 1/(i + j - 1) i=1,2,...N, j=1,2,..N
	//Due to the first element in a array being indexed at 0, we must add 1 to i and j in the equation to get to their correct values
	//this forces the equation to be A[i][j] = 1/(i + j + 1) i = 0...N-1, j = 0...N-1
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++) {
			hilbert[i][j] = 1.0L/(i + j + 1);
		}
	}

	return hilbert;
}

/**
 * calculateBinomial
 * Purpose: Calculate the binomial of choose k from n. This is used by the calculation of the Bn matrix
 *          This happens to use dynamic programming for the calculation, so that we don't neccisarily have to worry about overflow
 * Arguments: int n: the total number of objects to choose from
 *            int k: the number of objects to choose
 * Return value: the calculated binomial
 * See: http://www.csl.mtu.edu/cs4321/www/Lectures/Lecture%2015%20-%20Dynamic%20Programming%20Binomial%20Coefficients.htm
 *      (we had covered this way of calculating the binomial coeeficients in our CSCI 406 class, but I didn't quite remember how it was done, I had to look it up
 * Complexity: Time:  O(N^2)
 *             Space: O(N^2) // this could be lessened to be K, if you only store the previous row as well as the row we're currently on.  
 **/
long calculateBinomial(int n, int k) {
	//binom is a NxK matrix
	long** binom = new long*[n+1];
	for(int i = 0; i < (n+1); i++){
		binom[i] = new long[k+1];
		for(int j = 0; j <= k; j++){
			binom[i][j] = 0L;
		}
	}

	//binom[i][j] = binom[i-1][j-1] + binom[i-1][j]
	for(int i = 0; i <= n; i++) {
		for(int j = 0; j <= min(i,k); j++){
			if(i == j || j == 0){
				binom[i][j] = 1;
			}
			else { 
				binom[i][j] = binom[i-1][j-1] + binom[i-1][j];
			}
		}
	}
	return binom[n][k];
}

/**
 * calculateBtilde
 * Purpose: This is to calculate a value known as "B" with a tilde on top of it.
 *          This is a unique coefficient as far as I know, as I cannot think of something with a specific name for it.
 * Arguments: int n: the n argument to the B tilde function
 *            int k: the k argument to the B tilde argument
 * Returns: the value that B tilde should take.
 * Complexity: Time:  O(1)
 *             Space: O(1)
 */
long double calculateBtilde(int n, int k){
	long double btilde = 0.0L;
	btilde = cos((long double)n)/(cos((long double)k) * cos((long double)n - (long double)k));
	return btilde;
}

/**
 * This is main, the main entry point for an applicaiton.
 *
 **/
void main(){

}
