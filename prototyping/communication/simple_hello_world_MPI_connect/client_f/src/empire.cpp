#include "empire.h"

EMPIRE::EMPIRE() : portfile("server.port") {
message = 10101;
MPI_Init( NULL, NULL );  
}

EMPIRE::~EMPIRE(){}

void EMPIRE::Init(){
if (portfile.is_open())
{
     getline (portfile,port_name);      

     MPI_Comm_connect( (char*)(port_name.c_str()), MPI_INFO_NULL, 0, MPI_COMM_WORLD, &server );  
}
    else cout << "Unable to open file"; 
}

void EMPIRE::Send(){
      MPI_Send( &message, 1, MPI_INTEGER, 0, 0, server );     
}
void EMPIRE::Receive(){
      MPI_Recv( &message, 1, MPI_INTEGER, MPI_ANY_SOURCE, MPI_ANY_TAG, server, &status );
}
void EMPIRE::Disconnect(){
      cout << "RECEIVED:"<<message<<endl;
      MPI_Comm_size(MPI_COMM_WORLD, &size); 
      cout << "COMM size MPI_COMM_WORLD:"<<size<<endl;
      MPI_Comm_size(server, &size); 
      cout << "COMM size server:"<<size<<endl;
      MPI_Comm_test_inter(server, &isInterComm);
      cout << "IsInterComm client:"<<isInterComm<<endl;
      MPI_Comm_remote_size(server, &size); 
      cout << "COMM size remote server:"<<size<<endl;
      MPI_Comm_disconnect(&server);
      MPI_Finalize();
      portfile.close();
}