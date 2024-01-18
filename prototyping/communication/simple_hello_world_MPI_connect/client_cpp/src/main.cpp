#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#define MAX_DATA 1000
#define NUM_CLIENTS 5
int main( int argc, char **argv )
{
    MPI_Comm server;

    char port_name[MPI_MAX_PORT_NAME];
    int    size;
    MPI_Status status;
    int message, i, isInterComm;

    MPI_Init(0,0);
    fstream portfile;
    portfile.open ("server.port");
    //MPI_Init( &argc, &argv );

    string portName;
    portfile >> portName;
	cout << endl << endl << portName << endl << endl;

    int dummy =1;

    MPI_Comm_connect((char*) (portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server);
    MPI_Send(&dummy, 1, MPI_INTEGER, 0, 0, server);
    MPI_Recv(&dummy, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status);
    MPI_Comm_disconnect(&server);
}
