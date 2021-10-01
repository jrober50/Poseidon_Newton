   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE IO_Print_Results_Module                                               !##!
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
            ONLY :  idp, pi


USE Global_Variables_And_Parameters,&
            ONLY :  R_INNER,                    &
                    R_OUTER,                    &
                    Analytic_Solution

USE Poseidon_Main_Module, &
            ONLY :  Calculate_Potential_At_Location


IMPLICIT NONE


CONTAINS



!+101+##########################################################################!
!                                                                               !
!                                                     				!
!                                                                               !
!###############################################################################!
SUBROUTINE Print_Results()

INTEGER                                 ::  i
INTEGER                                 ::  NUM_SAMPLES
REAL(KIND = idp)                        ::  deltar, r, theta, phi
REAL(KIND = idp)                        ::  TMP, potential

110 FORMAT (11X,A1,16X,A17,10X,A12,10X,A19)
111 FORMAT (4X,A16,7X,A19,7X,A16,7X,A21)
112 FORMAT (ES22.15,3X,ES22.15,3X,ES22.15,3X,ES22.15)

NUM_SAMPLES = 20


deltar = (R_OUTER - R_INNER)/REAL(NUM_SAMPLES, KIND = idp)

theta = pi/2.0_idp
phi = pi/2.0_idp




WRITE(*,'(A,F4.2,A,F4.2,A)')"Results taken along ray, theta = ",theta/pi," Pi Radians, Phi = ",phi/pi," Pi Radians"
WRITE(*,110)"r","Analytic Solution","Code Results","Relative Difference"
WRITE(*,111)"----------------","-------------------","----------------","---------------------"

DO i = 0,NUM_SAMPLES

    r = i*deltar + R_INNER
    potential =  Calculate_Potential_At_Location( r )
    TMP = Analytic_Solution(r,theta,phi)


    IF ( ABS(TMP) == 0.0_idp) THEN

        WRITE(*,112)r,TMP,Potential,ABS(TMP-potential)/ABS(1e-16)

    ELSE

        WRITE(*,112)r,TMP,potential, ABS(TMP-potential)/ABS(TMP)

    END IF

END DO

PRINT*," "

END SUBROUTINE Print_Results








END MODULE  IO_Print_Results_Module
