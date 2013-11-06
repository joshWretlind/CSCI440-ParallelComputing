/***************************************
 * Author: Josh Wretlind
 * Date: 11/05/13
 * Purpose: This is the main file for my final project for CSCI 440
 *          (Performing SHA-3 hashing in parallel)
 ***************************************/


#include<iostream>
#include<string>
#include "mpi.h"
#include<cstdlib>
#include<bitset>

using namespace std;

int master = 0;
int totalSize;
int myRank;

/*******************************************************
 * convertStringTobits
 * Purpose: This method converts a string into straight binary bits
 * Arguments: string str: This is the string we wish to convert
 * Return value: bitset: the set of bits for the final result
 * Complexity:  Time: O()
 *             Space: O()
 * *****************************************************/
bitset<int> convertStringToBits(string str){
    
    int chunkPerWorker = ceil(((double)str.length())/((double)totalSize));
    int lowerBound = myRank*(chunkPerWorker);
    int upperBound = (rank+1)*(chunkPerWorker);
    if((rank+1) == totalSize){
	upperBound = str.length();
    }
    if(upperBound > str.length()){
	upperBound = str.length();
    }
    if(lowerBound > str.length()){
	lowerBound = str.length();
    }
    if(lowerBound == upperBound){
	return new bitset<0>();
    }
    bitset<8*(upperBound - lowerBound)> mainBitset;
    for(size_t i =lowerBound; i < upperBound; i++){
	bitset<8> currentChar = new bitset<8>(str.c_str()[i]);
	for(int j = 0; j < 8; j++){
	    mainBitset.set(8*i + j, currentChar[j]);
	} 
    }
    
    return mainBitset;
}

int main(int argc, char *argv[]){
    double startTime;
    double endTime;   
    
    //initialize Stuff
    MPI::Init(argc, argv);
    startTime = MPI::Wtime();
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    
    if(myRank == master){
    
    }
    
    int digest = atoi(argv[1]);
    string message = argv[2];
    
    cout << message << endl;
    
    
    //Finish things, clean up after ourselves.
    endTime = MPI::Wtime();
	MPI::Finalize();
}
