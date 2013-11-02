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
    double* rOfK = new double[size];
    
    for(int i = 0; i < size; i++){
        uniform_real_distribution<double> distribution(((double)myRank)/((double)i),i * myRank + 1.0);
        rOfK[i] = distribution(generator);
    }
    
    return rOfK;
 }


int main(int argc, char *argv[]){
    MPI::Init(argc, argv);
    
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    p = totalSize;
    cout << argv[1] << endl;
    
	MPI::Finalize();
}
