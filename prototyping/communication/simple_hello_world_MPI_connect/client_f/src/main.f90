!This is the client
!Using Intel MPI 4.0.3
!Starting the server with (mpdboot running) 
!mpiexec -np 1 ./server


PROGRAM main
! THIS TAKES CARE FOR EMPIRE API CROSS-COMPILATION CAPABILITIES
INTERFACE
    SUBROUTINE     EMPIRE_Init       BIND(C,NAME='EMPIRE_Init')
    END SUBROUTINE EMPIRE_Init

    SUBROUTINE     EMPIRE_Send       BIND(C,NAME='EMPIRE_Send')
    END SUBROUTINE EMPIRE_Send

    SUBROUTINE     EMPIRE_Receive    BIND(C,NAME='EMPIRE_Receive')
    END SUBROUTINE EMPIRE_Receive

    SUBROUTINE     EMPIRE_Disconnect BIND(C,NAME='EMPIRE_Disconnect')
    END SUBROUTINE EMPIRE_Disconnect
END INTERFACE
! HERE THE MAINPROGRAM STARTS
     print *, "Hello World from FORTRAN!"
     call EMPIRE_Init()
     call EMPIRE_Send()
     call EMPIRE_Receive()
     call EMPIRE_Disconnect()
END PROGRAM main




 
     

      
      

      
      















