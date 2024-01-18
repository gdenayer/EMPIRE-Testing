// This is the server
// Using Intel MPI 4.0.3
// Starting the server with (mpdboot running) 
//   mpiexec -np 1 ./server
   
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

#define MAX_DATA 1000
#define NUM_CLIENTS 5
int main( int argc, char **argv ) 
{ 
    MPI_Comm * InterClients = new MPI_Comm[NUM_CLIENTS];
    
    char port_name[MPI_MAX_PORT_NAME]; 
    int    size; 
    MPI_Status status; 
    int message, i, isInterComm; 
    
    ofstream portfile;
    portfile.open ("server.port");
    MPI_Init( &argc, &argv ); 
         
    MPI_Open_port(MPI_INFO_NULL, port_name); 
    cout <<"Server available at "<<port_name<<endl; 
    portfile << port_name << endl;
    portfile.close();
    
    
    for(i=0;i<NUM_CLIENTS;i++){
	MPI_Comm_accept( port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD, &InterClients[i] ); 
	cout << "Client " <<i<<" connected"<<endl;
	MPI_Recv( &message, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, InterClients[i], &status );
	cout << "MESSAGE: " << message <<endl;
	message = i+1;
	MPI_Send( &message, 1, MPI_INTEGER, 0, 0, InterClients[i] );
	cout << "SENDING..."<<endl;
	MPI_Comm_size(MPI_COMM_WORLD, &size); 
	cout << "COMM size MPI_COMM_WORLD:"<<size<<endl;
	MPI_Comm_size(InterClients[i], &size); 
	cout << "COMM size client "<<i<<":"<<size<<endl;
	MPI_Comm_test_inter(InterClients[i], &isInterComm);
	cout << "IsInterComm client "<<i<<":"<<isInterComm<<endl;
	MPI_Comm_remote_size(InterClients[i], &size); 
	cout << "COMM size remote client "<<i<<":"<<size<<endl;
	MPI_Comm_disconnect(&InterClients[i]);    
    }

    MPI_Finalize();

}
