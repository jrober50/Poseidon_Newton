   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Print_Setup_Module                                                 !##!
!##!                                                                         !##!
!##!_________________________________________________________________________!##!
!##!                                                                         !##!
!##!                                                                         !##!
!##!=========================================================================!##!
!##!                                                                         !##!
!##!    Contains:                                                            !##!
!##!                                                                         !##!
!##!                                                                         !##!
!###############################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !##########################################################################!


!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE constants, &
            ONLY :  idp


USE Global_Variables_And_Parameters,&
            ONLY :  Num_R_Elements,     &
                    Num_T_Elements,     &
                    Num_P_Elements,     &
                    Source_Degrees,     &
                    Degree,             &
                    L_Limit

IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                                     !
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Setup()


1400 FORMAT(/)
1401 FORMAT('------------- POSEIDON PARAMETERS --------------')

1402 FORMAT(/'------------- Expansion Parameters ------------')
1403 FORMAT('      Finite Element Degree = ',I3.1)
1404 FORMAT(' Spherical Harmonic L Limit = ',I3.1)

1405 FORMAT(/'------------- Mesh Parameters -----------------')
1406 FORMAT(' # Radial Elements = ',I5.1)
1407 FORMAT(' # Theta Elements  = ',I5.1)
1408 FORMAT(' # Phi Elements    = ',I5.1)


1410 FORMAT(/'------------- Quadrature Parameters -----------')
1411 FORMAT(' # Radial Quad Points = ',I3.1)
1412 FORMAT(' # Theta Quad Points  = ',I3.1)
1413 FORMAT(' # Phi Quad Points    = ',I3.1)


WRITE(*,1400)
WRITE(*,1401)
WRITE(*,1402)
WRITE(*,1403)Degree
WRITE(*,1404)L_LIMIT
WRITE(*,1405)
WRITE(*,1406)Num_R_Elements
WRITE(*,1407)Num_T_Elements
WRITE(*,1408)Num_P_Elements
WRITE(*,1410)
WRITE(*,1411)Source_Degrees(1)
WRITE(*,1412)Source_Degrees(2)
WRITE(*,1413)Source_Degrees(3)
WRITE(*,1400)






END SUBROUTINE Print_Setup








END MODULE  IO_Print_Setup_Module

