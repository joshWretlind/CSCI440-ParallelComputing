#include<iostream>
#include<cmath>
#include<vector>
#include<float.h>
#include<ctime>
#include<fstream>
#include<sstream>

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

long double factorial(int n){
	long double fact = 1;
	for(int i = 1; i <= n; i++){
		fact *= i;
	}
	return fact;
}
/**
 * calculateBinomial
 * Purpose: Calculates the binomial coefficient for choose k from n.
 * Arguments: int n: the total number of objects to choose from
 *            int k: the number of objects to choose
 * Return value: the calculated binomial
 * Complexity: Time:  O(N)
 **/
long double calculateBinomial(int n, int k) {
	return ((factorial(n))/(factorial(k) * factorial(n-k)));
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
 * Complexity: Time:  O(N)
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
 * Complexity: Time:  O(N)
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
 * Complexity: Time:  O(N)
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
 * Complexity: Time:  O(N^3)
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
 * Complexity: Time:  O(N^3)
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
		cout << "----------------------------" << endl;
		cout << "N: " << size << " Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << " Cn[N][N]: " << matrixC[size-1][size-1] << " T1Q1: " << difftime(endTime,startTime) << endl;
	}
}

/**
 * answerQuestion1Small
 * Purpose: This method is to provide all of the information needed to answer question 1 for smaller matrix sizes(1,2,4,8,16,32)
 * Complexity: Time:  O(N^3)
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

		cout << "----------------------------" << endl;
		cout << "N: " << size << " Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << " Cn[N][N]: " << matrixC[size-1][size-1] << " T1Q1: " << difftime(endTime,startTime) << endl;
	}
}

/**
 * answerQuestion12
 * Purpose: This method answers both parts of question 1 and 2 for the report
 * Complexity: Time:  O(N^3)
 *             Space: O(1)
 **/
void answerQuestion12() {
	cout << "------------------------------Answering questions 1 and 2-----------------------" << endl;
	answerQuestion1Small();
	answerQuestion2Large();
}

/**
 * matrixToString
 * Purpose: This method basically converts a matrix of long doubles to a string which can be written to disk.
 * Arguments: vector< vector<long double> > matrix: This is the matrix which you want to convert to a string
 * Return: The matrix converted to a string.
 * Complexity: Time:  O(N^2)
 *             Space: O(N^2)
 **/
string matrixToString(vector< vector<long double> > matrix){
	stringstream output;

	for(int i = 0; i < matrix.size();i++) {
		for(int j = 0; j < matrix.size(); j++) {
			output << matrix[i][j] << " ";
		}
		output << endl;
	}

	return output.str();
}
/**
 * writeMatriciesToDisk
 * Purpose: This creates and writes the A and B matricies to disk for 100, 200, 400, 800, 1600, 3200 sizes
 * Complexity: Time:  O(N^3)
 *             SpacE: O(N^2)
 **/
void writeMatriciesToDisk() {
	time_t startTime;
	time_t endTime;
	int size_base = 100;
	for(int i = 0; i < 6; i++){
		int size = size_base*pow(2.0, i);
		time(&startTime);
		vector< vector<long double> > matrixA = createHilbertMatrix(size);
		vector< vector<long double> > matrixB = calculateBinomialCoefficientMatrix(size);

		//create files for writing
		ofstream matrixAFile;
		ofstream matrixBFile;

		//make the file names
		stringstream matrixAFileName;
		matrixAFileName << "matrixA" << size << ".txt";

		stringstream matrixBFileName;
		matrixBFileName << "matrixB" << size << ".txt";
		
		//open files
		matrixAFile.open(matrixAFileName.str());
		matrixBFile.open(matrixBFileName.str());
		
		//write to the files
		matrixAFile << matrixToString(matrixA);
		matrixBFile << matrixToString(matrixB);

		matrixAFile.close();
		matrixBFile.close();

		time(&endTime);
		cout << "----------------------------" << endl;
		cout << "N: " << size << " T1Q1: " << difftime(endTime,startTime) << endl;
	}
}

/**
 * convertRowToMatrix
 * Purpose: This method is to convert a row from a string to a matrix.
 * Arguments: int i: The row to be converted/the values should be put into
 *            int n: the size of the matrix we are converting to
 *            vector< vector<long double> > &matrix: this is the matrix which we want to put the values into
 *            string input: the string of values we want to convert and put into the matrix
 * Complexity: Time:  O(N)
 *             Space: O(1)
 **/
void convertRowToMatrix(int i, int n, vector< vector<long double> > &matrix, string input){
	istringstream stringStream(input);
	string number = " ";
	for(int j = 0; getline(stringStream, number, ' '); j++){		
		stringstream  numberStream(number);
		numberStream >> matrix[i][j];
	}
}

/**
 * readMatriciesFromDiskAndMultiply
 * Purpose: This method reads in A and B matricies from the disk and then calculates C and then outputs some data
 * Complexity: Time:  O(N^3)
 *             Space: O(N^2)
 **/
void readMatriciesFromDiskAndMultiply(){
	time_t startTime;
	time_t endTime;
	int size_base = 100;
	for(int i = 0; i < 6; i++){
		int size = size_base*pow(2.0, i);
		time(&startTime);
		vector<long double> matrixARow(size, 0.0L);
		vector< vector<long double> > matrixA(size,matrixARow);

		vector<long double> matrixBRow(size, 0.0L);
		vector< vector<long double> > matrixB(size,matrixBRow);

		//Open the files
		ifstream matrixAFile;
		ifstream matrixBFile;
		string stringRow = " ";

		stringstream matrixAFileName;
		matrixAFileName << "matrixA" << size << ".txt";

		stringstream matrixBFileName;
		matrixBFileName << "matrixB" << size << ".txt";

		matrixAFile.open(matrixAFileName.str());
		matrixBFile.open(matrixBFileName.str());

		//Create the matricies from the files
		if(matrixAFile.is_open()){
			for(int i = 0; getline(matrixAFile,stringRow); i++){
				if(i != size){
					convertRowToMatrix(i,size,matrixA,stringRow);
				}
			}
		}
		//we don't need this file to stay open
		matrixAFile.close();
		
		if(matrixBFile.is_open()){
			for(int i = 0; getline(matrixBFile,stringRow); i++){
				if(i != size){
					convertRowToMatrix(i,size,matrixB,stringRow);
				}
			}
		}
		//close the file since it's no longer needed
		matrixBFile.close();

		vector< vector<long double> > matrixC = createMultipliedMatrix(matrixA, matrixB, size);

		time(&endTime);
		cout << "----------------------------" << endl;
		cout << "N: " << size << " Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << " Cn[N][N]: " << matrixC[size-1][size-1] << " T1Q1: " << difftime(endTime,startTime) << endl;
	}
}

/**
 * answerQuestion34
 * Purpose: This method answers both parts of question 3 and 4 for the report
 * Complexity: Time:  O(N^3)
 *             Space: O(1)
 **/
void answerQuestion34(){
	cout << "---------------------Writing and reading the Matricies to disk(questions 3 and 4)-----------------------------" << endl;
	writeMatriciesToDisk();
	readMatriciesFromDiskAndMultiply();
}

/**
 * answerQuestion5
 * Purpose: This is to output data for question #5 of the report(running it on 6400x6400 matricies
 * Complexity: Time:  O(N^3)
 *             Space: O(1)
 **/
void answerQuestion5(){
	time_t startTime;
	time_t endTime;
	int size = 100*pow(2.0, 6);
	time(&startTime);
	vector< vector<long double> > matrixA = createHilbertMatrix(size);
	vector< vector<long double> > matrixB = calculateBinomialCoefficientMatrix(size);
	vector< vector<long double> > matrixC = createMultipliedMatrix(matrixA, matrixB, size);
	time(&endTime);
	cout << "----------------------------" << endl;
	cout << "N: " << size << " Cn[N/2][N/2]: " << matrixC[(size)/2-1][(size/2)-1] << " Cn[N][N]: " << matrixC[size-1][size-1] << " T1Q1: " << difftime(endTime,startTime) << endl;
}

/**
 * This is main, the main entry point for an applicaiton.
 *
 **/
int main(){
	answerQuestion12();
	answerQuestion34();
}
