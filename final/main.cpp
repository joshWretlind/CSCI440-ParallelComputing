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
#include<cmath>
#include<vector>
#include<stdbool.h>

using namespace std;

int master = 0;
int totalSize;
int myRank;
int chunkPerWorker;
int lowerBound;
int upperBound;

/*******************************************************
 * convertStringTobits
 * Purpose: This method converts a string into straight binary bits
 * Arguments: string str: This is the string we wish to convert
 * Return value: bitset: the set of bits for the final result
 * Complexity:  Time: O()
 *             Space: O()
 * *****************************************************/
bool** convertStringToBits(string str){
    
    chunkPerWorker = ceil(((double)str.length())/((double)totalSize));
    lowerBound = myRank*(chunkPerWorker);
    upperBound = (myRank+1)*(chunkPerWorker);
    if((myRank+1) == totalSize){
	upperBound = str.length();
    }
    if(upperBound > str.length()){
	upperBound = str.length();
    }
    if(lowerBound > str.length()){
	lowerBound = str.length();
    }
    if(lowerBound == upperBound){
	return NULL;
    }
    bool** mainBitset = new bool*[upperBound - lowerBound];
    for(int i = 0; i < (upperBound - lowerBound); i++){
	mainBitset[i] = new bool[8];
    }
    for(int i =lowerBound; i < upperBound; i++){
	bitset<8> currentChar(str.c_str()[i]);
	for(int j = 0; j < 8; j++){
	    mainBitset[i - lowerBound][j] = currentChar[7-j];
	} 
    }
    
    return mainBitset;
}

bool* convertAndBroadcastBits(string message){
    bool** myBits = convertStringToBits(message);
    bool* messageInBinary = new bool[8*message.size()];

    if(myRank != master){
	for(int i = 0; i < (upperBound - lowerBound); i++){
	    MPI::COMM_WORLD.Send(myBits[i], 8, MPI_CHAR, master, myRank*chunkPerWorker + i); 
	}
    } else {
	for(int i = 0; i < chunkPerWorker; i++){
	    totalBits[i] = myBits[i];
	}
	MPI::Status myStatus;
	for(int i = (upperBound - lowerBound); i < message.length(); i++){
	    MPI::COMM_WORLD.Recv(totalBits[i], 8, MPI_CHAR, floor(((double)i)/((double)(upperBound - lowerBound))), i, myStatus);
	}
	for(int i = 0; i < message.size(); i++){
	    for(int j = 0; j < 8; j++){
		messageInBinary[8*i + j] = messageBits[i][j];
	    }
	}
	for(int i = 0; i < message.size(); i++){
	    delete[] totalBits[i];
	}
	delete[] totalBits;
    }

    MPI::COMM_WORLD.Bcast(messageInBinary[i], 8*message.size(), MPI_CHAR, master);
    


    for(int i = 0; i < chunkPerWorker; i++){
	delete[] myBits[i];
    }
    delete[] myBits;

    
    
    return messageInBinary;
}

int main(int argc, char *argv[]){
    double startTime;
    double endTime;   
    
    //initialize Stuff
    MPI::Init(argc, argv);
    startTime = MPI::Wtime();
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    int digest = atoi(argv[1]);
    string message = argv[2];
    int messageSize = 8*message.size();
    
    bool* messageInBinary = convertAndBroadcastBits(message);
    
 

    
    //Finish things, clean up after ourselves.
    endTime = MPI::Wtime();
    MPI::Finalize();
}
