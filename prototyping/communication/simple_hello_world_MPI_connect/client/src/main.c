// This is the client
// Using Intel MPI 4.0.3
// Starting the server with (mpdboot running) 
//   mpiexec -np 1 ./server
// and the client 
//   mpiexec -np 1 ./client 
// ===> This combination works ==> The Server prints "Connected"
// However using mpirun as recommended in the doc does not work
// Why does it not work? IS there a need to start some service 
// in order to support MPI_Publish_name and MPI_Lookup_name ?

#include <mpi.h> 
#include <stdio.h>
#include <string.h>
#define MAX_DATA 1000 
int main( int argc, char **argv ) 
{ 
    MPI_Comm server; 
    double buf[MAX_DATA]; 
    char port_name[MPI_MAX_PORT_NAME];
    char serv_name[256];
    strcpy( serv_name, "MyTest" );
    
    MPI_Init( &argc, &argv );    
    MPI_Lookup_name(serv_name,MPI_INFO_NULL,port_name);
    MPI_Comm_connect( port_name, MPI_INFO_NULL, 0, MPI_COMM_WORLD,  
                      &server ); 
    MPI_Comm_disconnect(&server);
    MPI_Finalize();
    return 0; 
}  
