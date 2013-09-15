#include<iostream>
#include<cmath>
#include<vector>
#include<float.h>
#include<ctime>

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
vector< vector<long double> > createHilbertMatrix(int n){
	//create a dynamic 2D array of N size
	vector<long double> row(n, 0.0L);
	vector< vector<long double> > hilbert(n,row);

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
long double calculateBinomial(int n, int k) {
	//binom is a NxK matrix
	vector<long double> row(k+1, 0.0L);
	vector< vector<long double> > binom(n+1,row);

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
 * calculateBtildeUniqueCoefficient
 * Purpose: This is to calculate a value that is used by a entry known as "B" with a tilde on top of it.
 *          This is a unique coefficient as far as I know, as I cannot think of something with a specific name for it.
 * Arguments: int n: the n argument to the B tilde function
 *            int k: the k argument to the B tilde argument
 * Returns: the value that B tilde should take.
 * Complexity: Time:  O(1)
 *             Space: O(1)
 */
long double calculateBtildeUniqueCoefficient(int n, int k){
	long double btilde = 0.0L;
	btilde = ((long double)(cos((long double)n)))/((long double)(cos((long double)k) * cos((long double)n - (long double)k)));
	return btilde;
}

/**
 * calculateBNormal
 * Purpose: this has the full computation/equation for the Bnormal(b entry without a tilde over it) value
 * Arguments: int n: This is the size of the matrix that you're computing values for
 *            int i: this is the row position for the entry you're computing the value for
 *            int j: this is the column position for the entry you're computing the value for
 * Returns: The value for BNormal
 * Complexity: Time:  O(N^2)
 *             Space: O(1)
 **/
long double calculateBNormal(int n, int i, int j) {
	long double bNormal = pow(-1.0L, i + j );
	bNormal *= (i + j - 1);
	bNormal *= calculateBinomial(n + i - 1, n - j );
	bNormal *= calculateBinomial(n + j - 1, n - i );
	bNormal *= pow(calculateBinomial(i + j - 2, i -1 ), 2.0L);
	return bNormal;
}

/**
 * calculateBTilde
 * Purpose: this has the full computation/equation for the Bnormal(b entry with a tilde over it) value
 * Arguments: int n: This is the size of the matrix that you're computing values for
 *            int i: this is the row position for the entry you're computing the value for
 *            int j: this is the column position for the entry you're computing the value for
 * Returns: The value for calculateBTilde
 * Complexity: Time:  O(N^2)
 *             Space: O(1)
 **/
long double calculateBTilde(int n, int i, int j) {
	long double bNormal = pow(-1.0L, i + j );
	bNormal *= (i + j - 1);
	bNormal *= calculateBtildeUniqueCoefficient(n + i -1, n - j );
	bNormal *= calculateBtildeUniqueCoefficient(n + j -1, n - i );
	bNormal *= pow(calculateBtildeUniqueCoefficient(i + j - 2, i-1), 2.0L);
	return bNormal;
}

/**
 * calculateBentry
 * Purpose: This function calculates what the values for each entry in the BinomialCoefficientMatrix should be
 * Arguments: int n: This is the size of the matrix
 *            int i: This is the row index for the entry for which you're calculating the value for
 *            int j: this is the column index for the entry for which you're calculating the value for
 * Return: the value which should exsist for this entry
 * Complexity: Time:  O(N^2)
 *             Space: O(1)
 **/
long double calculateBentry(int n, int i, int j) {
	long double bNormal = calculateBNormal(n,i,j);
	if(abs(bNormal) > pow(10.0L,50.0L) || bNormal != bNormal){
		return calculateBTilde(n,i,j);
	}
	return bNormal;
}

/**
 * calculateBinomialCoefficientMatrix
 * Purpose: this function creates the entire binomial(B) matrix
 * Arguments: int n: the size of the matrix you want to create
 * Returns: the created matrix
 * Complexity: Time:  O(N^4)
 *             Space: O(N^2)
 **/
vector< vector<long double> > calculateBinomialCoefficientMatrix(int n){
	vector<long double> row(n, 0.0L);
	vector< vector<long double> > binomialMatrix(n,row);

	for(int i = 0; i < n; i++) { 
		for(int j = 0; j < n; j++) {
			binomialMatrix[i][j] = calculateBentry(n, i +1 , j+1);
		}
	}

	return binomialMatrix;
}

/**
 * multiplyMatricies
 * Purpose: This method is to handle the multiplicaiton of two matricies together(in the format of AB). Note that this only handles square matricies at the moment.
 * Arguments: long double** matrixA: This is A in AB
 *            long double** matrixB: this is B in AB
 *            int n: this is the size of A and B
 * Return value: The multiplied matrix.
 * Complexity: Time:  O(N^3)
 *             SpacE: O(N^2)
 **/
vector< vector<long double> > multiplyMatricies(vector< vector<long double> > matrixA, vector< vector<long double> > matrixB, int n){
	vector<long double> row(n, 0.0L);
	vector< vector<long double> > resultant(n,row);

	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			long double resultEntry = 0.0L;
			for(int k = 0; k < n; k++){
				resultEntry += matrixA[i][k]*matrixB[k][j];
			}

			resultant[i][j] = resultEntry;
		}
	}

	return resultant;
}

/**
 * createMultipliedMatrix
 * Purpose: this is to create the full multiplied matrix, includeing the case where the entries are greater than 10^300
 * Arguments: long double** matrixA: this is the first matrix to multiply together with B, in the equation AB
 *            long double** matrixB: this is the second matrix to multiply together with A, in the equation AB
 *            int n: this is the size of A and B as far as matricies go/
 * Return Value: finalized multiplied matrix.
 * Complexity: Time:  O(N^3)
 *             Space: O(N^2)
 **/
vector< vector<long double> > createMultipliedMatrix(vector< vector<long double> > matrixA, vector< vector<long double> > matrixB, int n) {
	vector< vector<long double> > multipliedMatrix = multiplyMatricies(matrixA, matrixB, n);
	for(int i = 0; i < n; i++) {
		for(int j = 0; j < n; j++) {
			if(abs(multipliedMatrix[i][j]) > pow(10.0L, 300.0L) || multipliedMatrix[i][j] != multipliedMatrix[i][j]){
				multipliedMatrix[i][j] = cos((long double)(i+1) + (long double)(j+1));
			}
		}
	}

	return multipliedMatrix;
}

/**
 * printMatrix
 * Purpose: what this method does is it prints out a matrix to stdout
 * Arguments: vector<vector<long double>> matrix: This is the matrix which you want to print out
 * Complexity: Time:  O(n^2);
 *             Space: O(1);
 **/
void printMatrix(vector< vector<long double> > matrix){
	for(unsigned int i = 0; i < matrix.size(); i++){
		for(unsigned int j = 0; j < matrix.size(); j++){
			cout << matrix[i][j] << " ";
		}
		cout << endl;
	}
}

/**
 * answerQuestion2Large
 * Purpose: This is to answer the second question of our report for the larger matricies(100,200,400,800,1600,3200)
 * Complexity: Time:  O(N^4)
 *             Space: O(N^2)
 **/
void answerQuestion2Large(){
	time_t startTime;
	time_t endTime;
	int size_base = 100;
	for(int i = 0; i < 6; i++){
		int size = size_base*pow(2.0, i);
		time(&startTime);
		vector< vector<long double> > matrixA = createHilbertMatrix(size);
		vector< vector<long double> > matrixB = calculateBinomialCoefficientMatrix(size);
		vector< vector<long double> > matrixC = createMultipliedMatrix(matrixA, matrixB, size);
		time(&endTime);
		cout << "N: " << size << "Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << "Cn[N][N]: " << matrixC[size-1][size-1] << "T1Q1: " << difftime(startTime,endTime) << endl;
	}
}

/**
 * answerQuestion1Small
 * Purpose: This method is to provide all of the information needed to answer question 1 for smaller matrix sizes(1,2,4,8,16,32)
 * Complexity: Time:  O(N^4)
 *             Space: O(N^2)
 **/
void answerQuestion1Small() {
	time_t startTime;
	time_t endTime;
	for(int i = 1; i < 6; i++){
		int size = pow(2.0, i);

		time(&startTime);
		vector< vector<long double> > matrixA = createHilbertMatrix(size);
		vector< vector<long double> > matrixB = calculateBinomialCoefficientMatrix(size);
		vector< vector<long double> > matrixC = createMultipliedMatrix(matrixA, matrixB, size);
		time(&endTime);

		cout << "---------MatrixA-----------" << endl << endl;
		printMatrix(matrixA);
		cout << "---------MatrixB-----------" << endl << endl;
		printMatrix(matrixB);
		cout << "---------MatrixC-----------" << endl << endl;
		printMatrix(matrixC);

		cout << "N: " << size << " Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << " Cn[N][N]: " << matrixC[size-1][size-1] << " T1Q1: " << difftime(startTime,endTime) << endl;
	}
}

/**
 * answerQuestion2
 * Purpose: This method answers both parts of question 1 and 2 for the report
 * Complexity: Time:  O(N^4)
 *             Space: O(1)
 **/
void answerQuestion12() {
	answerQuestion1Small();
	answerQuestion2Large();
}
/**
 * This is main, the main entry point for an applicaiton.
 *
 **/
int main(){
	answerQuestion1();
}
