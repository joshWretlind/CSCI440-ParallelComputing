#include "mpi.h"
#include <iostream>
#include <vector>
#include <random>

using namespace std;

int master = 0;
int totalSize;
int myRank;

//Defined in the project description
int j;

//Total number of processors
int p;

/*****************************************************
 * generateRandomWeightedVector
 * Purpose: This method is to generate the size-dimentional random weighted
 *          vector(this is called r_k in the assignment)
 * Arguments: int size: this is the size that the vector should take
 * Returns: double*: This is the generated vector;
 * 
 * ***************************************************/
 double* generateRandomWeightedVector(int size){
    default_random_engine generator;
    time_t currentTime;
    time(&currentTime);
    currentTime += 100*myRank;
    generator.seed(currentTime);
    double* rOfK = new double[size];
    
    for(int i = 1; i <= size; i++){
        uniform_real_distribution<double> distribution(((double)myRank + 1)/((double)i),i * (myRank + 1) + 1.0);
        rOfK[i-1] = distribution(generator);
    }
    
    return rOfK;
 }


int main(int argc, char *argv[]){
    time_t startTime;
	time_t endTime;   
    //initialize Stuff
    MPI::Init(argc, argv);
    time(&startTime);
    	
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    p = totalSize;
    j = atoi(argv[1]);
    double** wMatrix = new double*[p];
    for(int i = 0; i < p; i++){
        wMatrix[i] = new double[j];
    }
    for(int i = 0; i < p; i++){
        for(int k = 0; k < j; k++){
            wMatrix[i][j] = 0;
        }
    }

    //generate random numbers
    double* rOfK = generateRandomWeightedVector(j);
    if(myRank != master){
        //If we aren't the master, send out values to master
        MPI::COMM_WORLD.Send(rOfK,j,MPI_DOUBLE,master,myRank);
    }
    else{
        //add our values to wMatrix
        wMatrix[0] = rOfK;
        MPI::Status my_status;
        //recieve values from all of the other hosts, add them to wMatrix
        for(int i = 1; i < totalSize; i++){
            MPI::COMM_WORLD.Recv(wMatrix[i],j,MPI_DOUBLE,i,i,my_status);
        }
    }
    
    double* normalizedVector = new double[p*j];
    if(myRank == master){
        //find the total/S
        double totalOfWeights;
        for(int i = 0; i < p; i++){
            for(int k = 0; k < j; k++){
                totalOfWeights += wMatrix[i][k];
            }
        }
        //create w/the normalized weight vector
        for(int i = 0; i < p; i++){
            for(int k = 0; k < j; k++){
                normalizedVector[i*j + k] = wMatrix[i][k]/totalOfWeights;
            }
        }
        
        //Sanity check on the normalized vector to make sure that it's valid
        double test = 0;
        for(int i = 0; i < (p*j); i++){
            test += normalizedVector[i];
        }
        //if it's outside of .1% of what we'd expect, abort
        if(abs(1.001 - test ) > 0.001){
            cout << "Aborting, normalized vector is too far off from what it should be(more than .1%)" << endl;
            MPI::COMM_WORLD.Abort(1);
        }
        
    }
    double* weightedVector = new double[j];
    MPI::COMM_WORLD.Scatter(normalizedVector, j, MPI_DOUBLE, weightedVector, j, MPI_DOUBLE, master); 
    
    cout << "My Rank " << myRank << " " << weightedVector[0] << endl;

    time(&endTime);
	MPI::Finalize();
}
