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


    double* rOfK = generateRandomWeightedVector(j);
    double firstValOfR = rOfK[0];
    if(myRank != master){
        MPI::COMM_WORLD.Send(rOfK,j,MPI_DOUBLE,master,myRank);
    }
    else{
        wMatrix[0] = rOfK;
        MPI::Status my_status;
        for(int i = 1; i < totalSize; i++){
            MPI::COMM_WORLD.Recv(wMatrix[i],j,MPI_DOUBLE,i,i,my_status);
        }
    }
    
    if(myRank == master){
        double totalOfWeights;
        for(int i = 0; i < p; i++){
            for(int k = 0; k < j; k++){
                totalOfWeights += wMatrix[i][k];
            }
        }
        double* normalizedVector = new double[p*j];
        for(int i = 0; i < p; i++){
            for(int k = 0; k < j; k++){
                normalizedVector[i*j + k] = wMatrix[i][k]/totalOfWeights;
            }
        }
        
        double test = 0;
        for(int i = 0; i < (p*j); i++){
            test += normalizedVector[i];
        }
        if(abs(test - 1.001) > 0.001){
            MPI::COMM_WORLD.Abort(1);
        }
    }
    time(&endTime);
	MPI::Finalize();
}
