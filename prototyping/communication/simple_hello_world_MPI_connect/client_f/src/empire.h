//EMPIRE_API.h
#include <mpi.h>
#include <iostream>
#include <fstream>
#include <string>
using namespace std;

class EMPIRE
{
  public:
     EMPIRE();                      
    ~EMPIRE();                     
     void Init(void);
     void Send(void);
     void Receive(void);
     void Disconnect(void);

 
  private:  
     MPI_Comm server;
     ifstream portfile;
     string port_name;
     int    size,isInterComm,message; 
     MPI_Status status; 
     
};
