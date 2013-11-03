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

double** generateXMatrix(int j, int p){
    double** xMat = new double*[j];
    for(int i = 0; i < j; i++){
        xMat[i] = new double[j*p];
    }
    
    default_random_engine generator;
    time_t currentTime;
    time(&currentTime);
    currentTime += 100*myRank;
    generator.seed(currentTime);
    
    for(int i = 1; i <= j; i++){
        for(int k = 1; k <= j*p; k++){
            uniform_real_distribution<double> distribution(-1.0*i*k,((double)i)/((double)k));
            xMat[i-1][k-1] = distribution(generator);
        }
    }
    
    return xMat;
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
        if(test > 1.001 || test < 0.999){
            cout << "Aborting, normalized vector is too far off from what it should be(more than .1%)" << endl;
            MPI::COMM_WORLD.Abort(1);
        }
        
    }
    
    MPI::COMM_WORLD.Bcast(normalizedVector, j*p, MPI_DOUBLE, master);
    
    double** xMatrix = generateXMatrix(j,p);
    
    double** cMatrix;
    if(myRank == master && j == 2 && p == 4){
        cMatrix = new double*[p*j];
        for(int i = 0; i < p*j; i++){
            cMatrix[i] = new double[p*j];
        }
    } else {
        cMatrix = new double*[j];
        for(int i = 0; i < j; i++){
            cMatrix[i] = new double[p*j];
        }
    }
    
    //Calculate all of the sample means for the x matrix we've generated
    //Store it in it's own vector for easy access.
    double* sampleMean = new double[j];
    for(int i = 0; i < j; i++){
        double sampleSum = 0;
        for(int k = 0; k < j*p; k++){
            sampleSum += xMatrix[i][k];
        }
        sampleMean[i] = sampleSum/((double)j*p);
    }
    
    double weightSum = 0;
    for(int i = 0; i < p*j; i++){
        weightSum += normalizedVector[i]*normalizedVector[i];
    }
    
    double** yMatrix = new double*[p*j];
    for(int i = 0; i < p*j; i++){
        yMatrix[i] = new double[p*j];
    }
    
    //calculate yMatrix
    for(int i = 0; i < j; i++){
        for(int k = 0; k < p*j; k++){
            yMatrix[i][k] = (normalizedVector[k] * (xMatrix[i][k] - sampleMean[i])) / sqrt((1 - weightSum));
        }
    }
    
    if(myRank != master){
        for(int i = 0; i < j; i++){
            MPI::COMM_WORLD.Send(yMatrix[i],p*j,MPI_DOUBLE,master,myRank*j + i);            
        }
    }
    else{
        for(int i = j; i < p*j; i++){
            MPI::Status myStatus;
            MPI::COMM_WORLD.Recv(yMatrix[i],p*j,MPI_DOUBLE,floor(((double)i)/j),i,myStatus);
        }
    }
    
    for(int i = 0; i < p*j; i++){
        MPI::COMM_WORLD.Bcast(yMatrix[i], j*p, MPI_DOUBLE, master);
    }
    
    double** yTranspose = new double*[p*j];
    for(int i = 0; i < p*j; i++){
        yTranspose[i] = new double[p*j];
    }
    
    //Calculate the transpose
    for(int i = 0; i < p*j; i++){
        for(int k = 0; k < p*j; k++){
            yTranspose[i][k] = yMatrix[k][i];
        }
    }
    
    //calculate C
    for(int i = 0; i < j; i++){
        for(int k = 0; k < p*j; k++){
            cMatrix[i][k] = 0;
            for(int l = 0; l < p*j; l++){
                cMatrix[i][k] += yMatrix[myRank*j + i][l] * yTranspose[l][myRank*j + k];
            }
        }
    }
    
    if(j == 2 && p == 4){
        if(myRank != master){
            for(int i = 0; i < j; i++){
                MPI::COMM_WORLD.Send(cMatrix[i],p*j,MPI_DOUBLE,master,myRank*j + i);            
            }
        } else {
            for(int i = j; i < p*j; i++){
                MPI::Status myStatus;
                MPI::COMM_WORLD.Recv(cMatrix[i],p*j,MPI_DOUBLE,floor(((double)i)/j),i,myStatus);
                for(int k = 0; k < p*j; k++){
                    cout << cMatrix[i][k] << " ";
                }
                cout << endl;
            }
            cout << "c[7] " << cMatrix[0][0] << endl;
            for(int i = 0; i < j*p; i++){
                for(int k = 0; k  < j*p; k++){
                    cout << cMatrix[i][k] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }
    
    time(&endTime);
	MPI::Finalize();
}
