/*           This file has been prepared for Doxygen automatic documentation generation.          */
/**************************************************************************************************/
// Question 1:
// When you comment out sleep(1) (line 52) the following program will have a deadlock.
// As MPI_Comm_accept is a blocking call and the port is already opened when MPI_Comm_connect is called
// in the second OMP section why it results in a deadlock?
/**************************************************************************************************/
// Question 2:
// Comment again sleep(1) and comment out MPI_Send and MPI_Recv (line 59 and 70).
// Why this results in an internal ABORT (Assertion failed in file...)?
/**************************************************************************************************/
#include <mpi.h>
#include <omp.h>
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

int main(int argc, char **argv) {

    // MPI Init
    int providedThreadSupport;
    MPI_Init_thread(NULL, NULL, MPI_THREAD_MULTIPLE, &providedThreadSupport);
    if (MPI_THREAD_MULTIPLE != providedThreadSupport) {
        cout << "WARNING: Requested MPI thread support is not guaranteed." << endl;
    }

    //Open Port
    string portName;
    char portNameChar[MPI_MAX_PORT_NAME];
    MPI_Open_port(MPI_INFO_NULL, portNameChar);
    portName = portNameChar;
    cout << "PortName:" << portName << endl;

    // Variables
    MPI_Comm interClient;
    MPI_Comm dummy;

    // Start two threads
    omp_set_num_threads(2);
#pragma omp parallel shared(interClient,portName,dummy)
    {
        /// Use OpemMP section construct for function parallelism (nowait = no barrier is called after the section ends)
#pragma omp sections
        {
#pragma omp section
            {
                cout << "Open the listening thread by using openmp" << endl;

                /// Question 1
                // sleep(1);
                MPI_Comm_accept((char*) (portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_WORLD,
                        &interClient);

                MPI_Status status;
                int message = 0;
                /// Question 2
                MPI_Recv(&message, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, interClient, &status);
                cout << "Received message: " << message << endl;
            }
#pragma omp section
            {
                cout << "Open the master thread by using openmp" << endl;
                MPI_Comm_connect((char*) (portName.c_str()), MPI_INFO_NULL, 0, MPI_COMM_SELF,
                        &dummy);

                int message = 10;
                /// Question 2
                MPI_Send(&message, 1, MPI_INT, 0, 0, dummy);
                MPI_Comm_disconnect(&dummy);

            }
        } /// End of sections

    } /// End of parallel section
    cout << "Stop openmp threading" << endl;
    MPI_Close_port((char*) (portName.c_str()));
    MPI_Finalize();
    return (0);
}
