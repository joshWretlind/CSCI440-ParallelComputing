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
    double startTime;
	double endTime;   
    //initialize Stuff
    MPI::Init(argc, argv);
    startTime = MPI::Wtime();
    
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    p = totalSize;
    j = atoi(argv[1]);
    double** wMatrix = new double*[p];
    for(int i = 0; i < p; i++){
        wMatrix[i] = new double[j];
    }

    //generate random numbers
    double* rOfK = new double[1*j];
    rOfK = generateRandomWeightedVector(j);
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
    
    for(int i = 0; i < p*j; i++){
        cout << normalizedVector[i] << " " ;
    }
    cout << endl;
    //calculate yMatrix
    for(int i = 0; i < j; i++){
        for(int k = 0; k < p*j; k++){
            yMatrix[i][k] = sqrt(normalizedVector[k])*((xMatrix[i][k] - sampleMean[i]))/(sqrt(1-weightSum));
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
                cMatrix[i][k] +=  yMatrix[myRank*j + i][l] * yTranspose[l][myRank*j + k];
            }
        }
    }
    
    double min = cMatrix[0][0];
    double max = cMatrix[0][0];
    
    //Format of this: [minimum, i_min, j_min, core]
    double minPair[4] = {min,0,0,((double)myRank)};
    double maxPair[4] = {max,0,0,((double)myRank)};
    
    for(int i =0; i < j; i++){
        for(int k = 0; k < p*j; k++){
            if(cMatrix[i][k] > max){
                max = cMatrix[i][k];
                maxPair[0] = max;
                maxPair[1] = myRank*j + i;
                maxPair[2] = k;
            }
            if(cMatrix[i][k] < min){
                min = cMatrix[i][k];
                minPair[0] = min;
                minPair[1] = myRank*j + i;
                minPair[2] = k;
            }
        }
    }
        
    //Container for the collected pairs
    
    double** collectedPairs;
    if(myRank == master){
        
    }
    if(myRank == master){
        cout << "FOR J = " << j << " AND P = " << p << "--------------------------------" << endl;
    }
    
    //Handle collection of the final C for the case of j=2 p=4
    if(j == 2 && p == 4){
        if(myRank != master){
            for(int i = 0; i < j; i++){
                MPI::COMM_WORLD.Send(cMatrix[i],p*j,MPI_DOUBLE,master,myRank*j + i);            
            }
        } else {
            for(int i = j; i < p*j; i++){
                MPI::Status myStatus;
                MPI::COMM_WORLD.Recv(cMatrix[i],p*j,MPI_DOUBLE,floor(((double)i)/j),i,myStatus);                
            }
            cout << "C Matrix ------------: " << endl;
            for(int i = 0; i < j*p; i++){
                for(int k = 0; k  < j*p; k++){
                    cout << cMatrix[i][k] << " ";
                }
                cout << endl;
            }
            cout << endl;
        }
    }
    
    if(myRank != master){
        MPI::COMM_WORLD.Send(&maxPair,4,MPI_DOUBLE,master,myRank);
    } else {
        double overallMax[4];
        
        double** collectedPairs = new double*[p];
        for(int i = 0; i < p; i++){
            collectedPairs[i] = new double[4];
        }
        
        overallMax[0] = maxPair[0];
        overallMax[1] = maxPair[1];
        overallMax[2] = maxPair[2];
        overallMax[3] = maxPair[3];
        
        collectedPairs[0][0] = maxPair[0];
        collectedPairs[0][1] = maxPair[1];
        collectedPairs[0][2] = maxPair[2];
        collectedPairs[0][3] = maxPair[3];

        for(int i = 1; i < p; i++){
            MPI::Status myStatus;
            MPI::COMM_WORLD.Recv(collectedPairs[i],4,MPI_DOUBLE,i,i,myStatus);
        }
        
        overallMax[0] = collectedPairs[0][0];
        
        for(int i = 0; i < p; i++){
            if(collectedPairs[i][0] > overallMax[0]){
                overallMax[0]= collectedPairs[i][0];
                overallMax[1]= collectedPairs[i][1];
                overallMax[2]= collectedPairs[i][2];
                overallMax[3]= collectedPairs[i][3];
            }
        }
        cout << "Cmax: " << overallMax[0] << " CoreMax: " << overallMax[3];
        cout << " i_max: " << overallMax[1] << " j_max " << overallMax[2] << endl;
        for(int i = 0; i < p; i++){
            delete[] collectedPairs[i];
        }
        delete[] collectedPairs;
    }
    
    if(myRank != master){
        MPI::COMM_WORLD.Send(&minPair,4,MPI_DOUBLE,master,myRank);
    } else {
        double** collectedPairs = new double*[p];
        
        for(int i = 0; i < p; i++){
            collectedPairs[i] = new double[4];
        }
        
        double overallMin[4];
        overallMin[0] = minPair[0];
        overallMin[1] = minPair[1];
        overallMin[2] = minPair[2];
        overallMin[3] = minPair[3];
        
        collectedPairs[0][0] = minPair[0];
        collectedPairs[0][1] = minPair[1];
        collectedPairs[0][2] = minPair[2];
        collectedPairs[0][3] = minPair[3];
        for(int i = 1; i < p; i++){
            MPI::Status myStatus;
            MPI::COMM_WORLD.Recv(collectedPairs[i],4,MPI_DOUBLE,i,i,myStatus);
        }
        
        overallMin[0] = collectedPairs[0][0];
        
        for(int i = 0; i < p; i++){
            if(collectedPairs[i][0] < overallMin[0]){
                overallMin[0] = collectedPairs[i][0];
                overallMin[1] = collectedPairs[i][1];
                overallMin[2] = collectedPairs[i][2];
                overallMin[3] = collectedPairs[i][3];
            }
        }
        cout << "Cmin: " << overallMin[0] << " CoreMin: " << overallMin[3];
        cout << " i_min: " << overallMin[1] << " j_min " << overallMin[2] << endl;
        
        for(int i = 0; i < p; i++){
            delete[] collectedPairs[i];
        }
        delete[] collectedPairs;
    }
    endTime = MPI::Wtime();
	MPI::Finalize();
    
    if(j != 2 || p != 4){
        cout << "My Rank was: "  << myRank << " and I took " << endTime-startTime << " seconds" << endl ; 
    }
    
    //clean up our memory
    delete[] rOfK;
    for(int i = 0; i < j; ++i){
        delete[] xMatrix[i];
    }
    delete[] xMatrix;
    for(int i = 0; i < p*j; i++){
        delete[] yMatrix[i];
    }
    delete[] yMatrix;
    for(int i = 0; i < p*j; i++){
        delete[] yTranspose[i];
    }
    delete[] yTranspose;
    delete[] normalizedVector;
    delete[] sampleMean;
}
