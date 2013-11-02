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
    currentTime + 100*myRank;
    generator.seed(currentTime);
    double* rOfK = new double[size];
    
    for(int i = 0; i < size; i++){
        uniform_real_distribution<double> distribution(((double)myRank)/((double)i),i * myRank + 1.0);
        rOfK[i] = distribution(generator);
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
    const int j_size = j;
    const int p_size = p;
    double** wMatrix = new double[p_size][j_size];


    double* rOfK = generateRandomWeightedVector(j);
    
    if(myRank != master){
        MPI::COMM_WORLD.send(rOfK,j,MPI_DOUBLE,master,myRank);
    }
    else{
        MPI::Status my_status;
        for(int i = 1; i < totalSize; i++){
            MPI::COMM_WORLD.Recv(wMatrix[i],j,MPI_DOUBLE,i,i,my_status);
        }
    }
    delete rOfK;
    
    if(myRank == master){
        for(int i = 0; i < p; i++){
            for(int k = 0; k < j; j++){
                cout << wMatrix[i][k] << " ";
            }
            cout << endl;
        }
    }
    time(&endTime);
	MPI::Finalize();
}
