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
//l is the power of the word length size. This directly impacts most of the
//RAM requirements of this code
int l = 6;
//w is the word size. w=2^l
int w;
//c is the word "capacity". this is twice the imput number. Ex: SHA3-512 has c=1024
int c;
//r is the rate per permutation we need to do. r = 25*w - c
int r;
//Number of Rounds per each block, 2*l + 12
int numOfRounds;
//b is the width of each state permutation, 25*w
int b;
//RC is an array of 64-bit round constants
long *rc;
int cycleOffset[5][5];

bool*** b;

int messageSize;
int paddedSize;
/*******************************************************
 * convertStringTobits
 * Purpose: This method converts a string into straight binary bits
 * Arguments: string str: This is the string we wish to convert
 * Return value: bitset: the set of bits for the final result
 * Complexity:  Time: O(N)
 *             Space: O(N)
 *******************************************************/
bool** convertStringToBits(string str){
    
    //Figure out how mand chunks per worker there should be
    //As well as figure out what the lower and upper bounds should be
    chunkPerWorker = ceil(((double)str.length())/((double)totalSize));
    lowerBound = myRank*(chunkPerWorker);
    upperBound = (myRank+1)*(chunkPerWorker);
    
    //Handle the cases on the edge of the string
    if((myRank+1) == totalSize){
	upperBound = str.length();
    }
    if(upperBound > str.length()){
	upperBound = str.length();
    }
    if(lowerBound > str.length()){
	lowerBound = str.length();
    }
    
    //If our lowerBound == upperBound, we shouldn't do any conversion
    //Let this worker short circuit
    if(lowerBound == upperBound){
	return NULL;
    }
    //Initiallize our bits
    bool** mainBitset = new bool*[upperBound - lowerBound];
    for(int i = 0; i < (upperBound - lowerBound); i++){
	mainBitset[i] = new bool[8];
    }
    //Convert the string into bits
    for(int i =lowerBound; i < upperBound; i++){
	bitset<8> currentChar(str.c_str()[i]);
	for(int j = 0; j < 8; j++){
	    mainBitset[i - lowerBound][j] = currentChar[j];
	} 
    }
    
    return mainBitset;
}

/***************************************************
 * convertAndBroadcastBits
 * Purpose: This function basically is the wrapper for converting the
 *          message into bits, padding it, and then sending back the padded
 *          message back to all the hosts.
 * Input: string message: This is the message comming in from the command line
 * Output: an array representing the padded output message
 * Complexity: Time:  O(N)
 *             Space: O(N)
 **************************************************/
bool* convertAndBroadcastBits(string message){
    int padding = ceil(((double)messageSize)/((double)r));
    paddedSize = messageSize*r;
    
    bool** myBits = convertStringToBits(message);
    bool* messageInBinary = new bool[paddedSize];
    
    if(myRank != master){
	//Send our stuff to master
	for(int i = 0; i < (upperBound - lowerBound); i++){
	    MPI::COMM_WORLD.Send(myBits[i], 8, MPI_CHAR, master, myRank*chunkPerWorker + i); 
	}
    } else {
	//Initialize our message array
	bool** totalBits = new bool*[message.size()];
	for(int i = 0; i < message.size(); i++){
	    totalBits[i] = new bool[8];
	}
	//Set the first rows to the bits we calculated
	for(int i = 0; i < chunkPerWorker; i++){
	    totalBits[i] = myBits[i];
	}
	//recieve all the bits from the hosts/slaves
	MPI::Status myStatus;
	for(int i = (upperBound - lowerBound); i < message.length(); i++){
	    MPI::COMM_WORLD.Recv(totalBits[i], 8, MPI_CHAR, floor(((double)i)/((double)(upperBound - lowerBound))), i, myStatus);
	}
	//Flatten our array from a 2D array to a 1D array
	for(int i = 0; i < message.size(); i++){
	    for(int j = 0; j < 8; j++){
		messageInBinary[8*i + j] = totalBits[i][j];
	    }
	}
	//Handle padding
	if((paddedSize - messageSize) == 1) {
	    messageInBinary[paddedSize -1 ] = true;
	} else if((paddedSize - messageSize) > 1) {
	    messageInBinary[messageSize] = true;
	    for(int i = messageSize + 1; i < (paddedSize - 1); i++) {
		messageInBinary[i] = false;
	    }
	    messageInBinary[paddedSize - 1] = true; 
	}
	//Clean up our memory
	for(int i = 0; i < message.size(); i++){
	    delete[] totalBits[i];
	}
	delete[] totalBits;
    }
    //every one recieves the final message
    MPI::COMM_WORLD.Bcast(messageInBinary, paddedSize, MPI_CHAR, master);

    return messageInBinary;
}

/*********************************
 * 
 * 
 * 
 *********************************/
bool* rotate(bool* bitsToRotate, int howMuchToRotate,int lengthOfBitsToRotate){
    bool* rotated = new bool[lengthOfBitsToRotate];
    for(int i = 0; i < lengthOfBitsToRotate; i++){
	rotated[i] = bitsToRotate[(i+howMuchToRotate)%lengthOfBitsToRotate];
    }
    return rotated;
}

/*********************************
 * 
 * 
 * 
 *********************************/
void thetaStep(int totalKeccakSize, bool*** tempState, long rc){
    bool** c = new bool*[5];
    bool** d = new bool*[5];
    for(int i = 0; i < 5; i++){
	c[i] = new bool[w];
	d[i] = new bool[w];
    }
    
    for(int i = 0; i < 5; i++){
	for(int k = 0; k < w; k++){
	    c[i][k] = tempState[i][0][k] ^ tempState[i][1][k] ^
		      tempState[i][2][k] ^ tempState[i][3][k] ^
		      tempState[i][4][k];
	}
    }
    
    for(int i = 0; i < 5; i++){
	for(int j = 0; j < w; j++){
	    int backPos = i - 1 ;
	    if(backPos == -1){
		backPos = 4;
	    }
	    d[i][j] = c[backPos][j] ^ c[(i + 1)%5][(j+1)%w];
	}
    }
    
    for(int i = 0; i < 5; i++){
	for(int j = 0; j < 5; j++){
	    for(int k = 0; k < w; k++){
		tempState[i][j][k] = tempState[i][j][k] ^ d[i][k];
	    }
	}
    }
    
    for(int i = 0; i < 5; i++){
	delete[] c[i];
	delete[] d[i];
    }
    delete[] c;
    delete[] d;
}

/**************************************
 * 
 * 
 * 
 * 
 **************************************/
void piAndRhoSteps(int totalKeccakSize, bool*** tempState, long rc){
    b = new bool**[5];
    for(int i = 0; i < 5; i++){
	b[i] = new bool*[5];
	for(int j = 0; j < 5; j++){
	    b[i][j] = new bool[w];
	}
    }
    
    for(int i = 0; i < 5; i++){
	for(int j = 0; j < 5; j++){
	    bool* rotated = rotate(tempState[i][j],cycleOffset[i][j],w);
	    for(int k = 0; k < w; k++){
		b[j][(2*i + 3*y)%5][k] = rotated[k];
	    }
	}
    }
}

/**************************************
 * 
 * 
 * 
 * 
 **************************************/
void chiStep(int totalKeccakSize, bool*** tempState, long rc){
    for(int i = 0; i < 5; i++){
	for(int j = 0; j < 5; j++){
	    for(int k = 0; k < w; k++){
		tempState[i][j][k] = b[i][j][k] ^ ((!b[(i+1)%5][j][k]) & b[(i+2)%5][j][k]);
	    }
	}
    }
}

/**************************************
 *
 * 
 * 
 * 
 **************************************/
void iotaStep(int totalKeccakSize, bool*** tempState, long rc){
    bitset<64> rcVal(rc);
    for(int i = 0; i < w; i++){
	tempState[0][0][i] = rcVal[i] ^ tempState[0][0][i];
    }
}

/**************************************
 * keccakRound
 * 
 * 
 * 
 **************************************/
void keccakRound(int totalKeccakSize, bool*** tempState, long rc){
    thetaStep    (totalKeccakSize, tempState, rc);
    piAndRhoSteps(totalKeccakSize, tempState, rc);
    chiStep      (totalKeccakSize, tempState, rc);
    iotaStep     (totalKeccakSize, tempState, rc);
}

/***************************************
 * keccakSponge
 * 
 * 
 **************************************/
void keccakSponge(int totalKeccakSize, bool*** absorbedState){
    for(int i = 0; i < numOfRounds; i++){
	keccakRound(b,absorbedState,rc[i]);
    }
}
/**************************************
 * 
 * 
 * 
 **************************************/
void setupRoundConstants(){
    rc = new long[24];
    
    rc[0] = 0x1;
    rc[1] = 0x8082;
    rc[2] = 0x800000000000808A;
    rc[3] = 0x8000000080008000;
    rc[4] = 0x808B;
    rc[5] = 0x80000001;
    rc[6] = 0x8000000080008081;
    rc[7] = 0x8000000000008009;
    rc[8] = 0x8A;
    rc[9] = 0x88;
    rc[10] = 0x80008009;
    rc[11] = 0x8000000A;
    rc[12] = 0x8000808B;
    rc[13] = 0x800000000000008B;
    rc[14] = 0x8000000000008089;
    rc[15] = 0x8000000000008003;
    rc[16] = 0x8000000000008003;
    rc[17] = 0x8000000000008003;
    rc[18] = 0x800A;
    rc[19] = 0x800000008000000A;
    rc[20] = 0x8000000080008081;
    rc[21] = 0x8000000000008080;
    rc[22] = 0x80000001;
    rc[23] = 0x8000000080008008;
}
/*********************************
 * 
 * 
 * 
 * 
 * *******************************/
void setupCyclicConstants(){
    int cycleOffset[5][5] = {{0,1,62,28,27},
			     {36,44,6, 55,20},
                             {3,13,43,25,39},
                             {41,45,15,21,8},
                             {18,2,61,56,14}};
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
    messageSize = 8*message.size();
    
    //Set up some constants for keccak
    c = 2*digest;
    w = pow(2.0,l);
    r = 25*w - c;
    b = 25*w;
    numOfRounds = 12 + 2*l;
    
    setupRoundConstants();
    setupCyclicConstants();
    
    bool* messageInBinary = convertAndBroadcastBits(message);
        
    //Finish things, clean up after ourselves.
    endTime = MPI::Wtime();
    MPI::Finalize();
}
