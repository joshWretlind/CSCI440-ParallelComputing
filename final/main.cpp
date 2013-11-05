#include<iostream>
#include "mpi.h"

using namespace std;

int master = 0;
int totalSize;
int myRank;

int main(int argc, char *argv[]){
	double startTime;
	double endTime;   
    
    //initialize Stuff
    MPI::Init(argc, argv);
    startTime = MPI::Wtime();
    totalSize = MPI::COMM_WORLD.Get_size();
    myRank = MPI::COMM_WORLD.Get_rank();
    
    
    //Finish things, clean up after ourselves.
    endTime = MPI::Wtime();
	MPI::Finalize();
}
