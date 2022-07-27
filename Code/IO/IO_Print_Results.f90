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

USE Global_Variables_And_Parameters, &
            ONLY :  DEGREE,                 &
                    L_LIMIT,                &
                    R_INNER,                &
                    R_OUTER,                &
                    NUM_R_ELEMENTS,         &
                    NUM_T_ELEMENTS,         &
                    NUM_P_ELEMENTS,         &
                    NUM_R_NODES,            &
                    rlocs,                  &
                    tlocs,                  &
                    plocs,                  &
                    Coefficient_Vector


USE Additional_Functions_Module, &
            ONLY :  Lagrange_Poly,                  &
                    Spherical_Harmonic,             &
                    Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LGL_Quadrature

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








!+300+##################################################################!
!                                                                       !
!       CALC_POTENTIAL - Use coefficents to calculate the potential.    !
!                                                                       !
!#######################################################################!
FUNCTION Calculate_Potential_At_Location( r_input, theta_input, phi_input)

REAL(KIND = idp)                                               ::  Calculate_Potential_At_Location


REAL(KIND = idp),               INTENT(IN)                      ::  r_input
REAL(KIND = idp),               INTENT(IN), OPTIONAL            ::  theta_input, phi_input




INTEGER                                 ::  l, m, re, rd
REAL(KIND = idp)                        ::  r, theta, phi
REAL(KIND = idp)                        ::  r_tmp
REAL(KIND = idp), DIMENSION(0:DEGREE)   ::  LagP
REAL(KIND = idp), DIMENSION(0:DEGREE)   ::  xlocP, weightP


COMPLEX(KIND = idp)                     ::  potential



r = r_input

if ( PRESENT(theta_input) ) THEN

    theta = theta_input

ELSE

    theta = 0.0_idp

END IF


IF (PRESENT(phi_input) ) THEN

    phi = phi_input

ELSE

    phi = 0.0_idp

END IF






potential = 0.0_idp

IF (r == rlocs(0)) THEN

    DO l = 0,L_LIMIT

        DO m = -l,l

            potential = potential + Coefficient_Vector(0,m,l)*Spherical_Harmonic(l,m,theta,phi)

        END DO

    END DO

ELSE

    CALL Initialize_LGL_Quadrature(DEGREE,xlocP,weightP)


    DO re = 0, NUM_R_ELEMENTS-1

        IF (r > rlocs(re) .AND. r <= rlocs(re+1)) THEN

            r_tmp = Map_To_X_Space(rlocs(re),rlocs(re+1),r)

            LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)


            DO l = 0,L_LIMIT

                DO m = -l,l

                    DO rd = 0,DEGREE


                        potential = potential + Coefficient_Vector(re*DEGREE + rd, m, l)*Spherical_Harmonic(l,m,theta,phi)*LagP(rd)


                    END DO

                END DO

            END DO


        END IF

    END DO



    IF (    r > rlocs(NUM_R_ELEMENTS)   ) THEN



        r_tmp = Map_To_X_Space(rlocs(NUM_R_ELEMENTS-1),rlocs(NUM_R_ELEMENTS),r)

        LagP = Lagrange_Poly(r_tmp,DEGREE,xlocP)


        DO l = 0,L_LIMIT

            DO m = -l,l

                DO rd = 0,DEGREE


                    potential = potential + Coefficient_Vector(NUM_R_NODES-1, m, l)    &
                                            * Spherical_Harmonic(l,m,theta,phi)                      &
                                            * LagP(rd)


                END DO

            END DO

        END DO



    END IF

END IF



Calculate_Potential_At_Location = REAL( potential )




END FUNCTION Calculate_Potential_At_Location








END MODULE  IO_Print_Results_Module
