   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
MODULE Error_Analysis_Functions                                                     !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains functions and subroutines to calculate L_ONE, L2, and L_Infinity   !##!
!##!    norm errors.                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!   +101+   L_ONE_ERROR                                                          !##!
!##!   +102+   L_TWO_ERROR_3D                                                       !##!
!##!   +103+   L_INF_ERROR                                                          !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!                                                                                !##!
!######################################################################################!
 !\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/!
   !################################################################################!

!*D*================================!
!                                   !
!           Dependencies            !
!                                   !
!===================================!
USE constants, &
                ONLY : idp, pi


USE Global_Variables_And_Parameters, &
                ONLY :  R_INNER, R_OUTER,                               &
                        NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                        L_LIMIT, DEGREE, Source_Function_Flag,          &
                        rlocs, tlocs, plocs,                            &
                        Analytic_Solution



USE Additional_Functions_Module, &
                ONLY :  Map_From_X_Space, Initialize_LG_Quadrature,            &
                        Initialize_LGL_Quadrature



USE Poseidon_Main_Module, &
                ONLY :  Calculate_Potential_At_Location


IMPLICIT NONE






!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS





 !+101+########################################################################!
!                                                                               !
!   L_ONE_ERROR - Calculates the L_one Error.  Good for analyzing error along   !
!                   a radial line outward.                                      !
!                                                                               !
 !#############################################################################!
SUBROUTINE L_ONE_ERROR(L_ONE_ERR,CUR, ITERS)

INTEGER, INTENT(IN)                                             :: CUR,ITERS
REAL(KIND = idp),DIMENSION(0:ITERS), INTENT(INOUT)              :: L_ONE_ERR
INTEGER                                                         :: ORDER, i


REAL(KIND = idp)                                       :: Diff_Int, theta, phi
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)            :: xlocs, rlocs, weights, Exp, Sol
REAL(KIND  = idp)                                    :: potential

!theta = 0.0_idp
!phi = 0.0_idp
theta = pi/2.0_idp
phi = pi/2.0_idp

ORDER = 30


ALLOCATE(xlocs(0:ORDER),rlocs(0:ORDER),weights(0:ORDER),Exp(0:ORDER),Sol(0:ORDER))


CALL Initialize_LGL_Quadrature(ORDER,xlocs,weights)

rlocs = Map_From_X_Space(R_INNER, R_OUTER, xlocs)

DO i = 0,ORDER

    potential = Calculate_Potential_At_Location( rlocs(i), theta, phi)


    Exp(i) = potential
    Sol(i) = Analytic_Solution(rlocs(i),theta,phi)

    !print*,Exp(i), Sol(i), Exp(i)-Sol(i)

END DO



Diff_Int = DOT_PRODUCT(weights,ABS(Exp-Sol))



DEALLOCATE(xlocs,rlocs,weights,Exp,Sol)


L_ONE_ERR(CUR) = Diff_Int/(R_OUTER-R_INNER)

END SUBROUTINE L_ONE_ERROR



 !+102+########################################################################!
!                                                                               !
!   L_TWO_ERROR_3D - Subroutine that calculates the L_Two erro integral over    !
!                       the entire computational domain.  Also reports maximum  !
!                       and minimum relative error from quadrature points       !
!                                                                               !
 !#############################################################################!
SUBROUTINE L_TWO_ERROR_3D(InnerTime)


REAL(KIND = idp), INTENT(IN)                        :: InnerTime

INTEGER, PARAMETER                                  :: OUTPUT_PARAMS = 1
INTEGER, PARAMETER                                  :: OUTPUT_TO_FILE = 1



CHARACTER(len = 55)                                 ::  outfile


INTEGER                                             ::  re, rd, &
                                                        te, td, &
                                                        pe, pd



INTEGER                                             ::  R_Degree, P_Degree, T_Degree

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  R_xlocs, R_locs, R_weights, &
                                                        P_xlocs, P_locs, P_weights, &
                                                        T_xlocs, T_locs, T_weights


REAL(KIND = idp)                                    ::  deltar, deltat, deltap

REAL(KIND = idp)                                    ::  R_PRE, T_PRE


REAL(KIND = idp)                                 ::  TMP_VALUE_EXP, TMP_VALUE_SOL

REAL(KIND = idp)                                    ::  TMP_INTEGRAL,                   &
                                                        TMP_DIFFERENCE,                 &
                                                        TMP_MAX,                        &
                                                        TMP_MIN

R_Degree = 10
P_Degree = 10
T_Degree = 10


ALLOCATE(R_xlocs(1:R_Degree), R_locs(1:R_Degree), R_weights(1:R_Degree))
ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))


CALL Initialize_LG_Quadrature(R_Degree, R_xlocs, R_weights)
CALL Initialize_LGL_Quadrature(P_Degree, P_xlocs, P_weights)
CALL Initialize_LGL_Quadrature(T_Degree, T_xlocs, T_weights)


TMP_INTEGRAL = 0.0_idp
TMP_MAX = TINY(0.0_idp)
TMP_MIN = HUGE(0.0_idp)



DO re = 0, NUM_R_ELEMENTS-1

    !-------------------------------------------------------!
    !                                                       !
    !   Map xspace quadratures points to real radial space  !
    !   Caclulate delta r                                   !
    !                                                       !
    !-------------------------------------------------------!
    R_locs = Map_From_X_Space(rlocs(re),rlocs(re+1), R_xlocs)
    deltar = (rlocs(re+1) - rlocs(re))



    DO te = 0,NUM_T_ELEMENTS-1


        !-------------------------------------------------------!
        !                                                       !
        !   Map xspace quadratures points to real theta space   !
        !   Caclulate delta t                                   !
        !                                                       !
        !-------------------------------------------------------!
        T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), T_xlocs)
        deltat = tlocs(te+1)-tlocs(te)


        DO pe = 0,NUM_P_ELEMENTS-1



            !---------------------------------------------------!
            !                                                   !
            !   Map xspace quadratures points to real phi space !
            !   Caclulate delta p                               !
            !                                                   !
            !---------------------------------------------------!
            P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), P_xlocs)
            deltap = plocs(pe+1) - plocs(pe)





            DO rd = 1,R_Degree


                !---------------------------------------------------!
                !                                                   !
                ! Precalculate radial portion of integral           !
                !                                                   !
                !---------------------------------------------------!
                R_PRE =  R_locs(rd)*R_locs(rd)*R_weights(rd)*deltar/2.0_idp





                DO td = 0,T_Degree


                    !---------------------------------------------------!
                    !                                                   !
                    ! Precalculate theta portion of integral            !
                    !                                                   !
                    !---------------------------------------------------!
                    T_PRE = sin(T_locs(td)) * T_weights(td) *  deltat/2.0_idp





                    DO pd = 0,P_Degree

                        !
                        ! Find potential given by method at location
                        !
                        TMP_VALUE_EXP = Calculate_Potential_At_Location(R_locs(rd),T_locs(td),P_locs(pd))


                        !
                        ! Find potential given by analytic solution at location
                        !
                        TMP_VALUE_SOL = Analytic_Solution(R_locs(rd),T_locs(td),P_locs(pd))


                        !
                        ! Find difference between method and analytic solutions
                        !
                        TMP_DIFFERENCE = ABS(TMP_VALUE_EXP - TMP_VALUE_SOL)

                        IF (TMP_VALUE_SOL .NE. 0.0_idp) THEN


                            TMP_DIFFERENCE = TMP_DIFFERENCE/TMP_VALUE_SOL

                        END IF


                        !---------------------------------------------------!
                        !                                                   !
                        ! Perform L2 Error Integral                         !
                        !                                                   !
                        !---------------------------------------------------!
                        TMP_INTEGRAL =  TMP_INTEGRAL                        &
                                        + sqrt(TMP_DIFFERENCE*TMP_DIFFERENCE)                    &
                                        * R_PRE                             &
                                        * T_PRE                             &
                                        * P_weights(pd)* deltap/2.0_idp





                        !---------------------------------------------------!
                        !                                                   !
                        ! Calculate Relative Error and store max/min value  !
                        !                                                   !
                        !---------------------------------------------------!

                        IF (TMP_DIFFERENCE > TMP_MAX) THEN

                            TMP_MAX = TMP_DIFFERENCE


                        ELSE IF (TMP_DIFFERENCE < TMP_MIN) THEN

                            TMP_MIN = TMP_DIFFERENCE

                        END IF





                    END DO ! pd loop

                END DO ! td loop

            END DO ! rd loop

        END DO ! pe loop

    END DO ! te loop

END DO ! re loop



IF (OUTPUT_TO_FILE == 1) THEN

    WRITE(outfile,'(A,i1,A,i1,A)'),'OUTPUT/NEWTON/ERROR/ERROR_OUTPUT_L',L_LIMIT,'_DEG',DEGREE,'.dat'
    OPEN(unit = 101, file = outfile )


    WRITE(101,*),"For a run with paramters: "
    WRITE(101,*),"-------------------------------------------------"
    WRITE(101,*),"Source Function               : ",Source_Function_Flag
    WRITE(101,*),"-------------------------------------------------"
    WRITE(101,*),"Number of Radial Elements     : ", NUM_R_ELEMENTS
    WRITE(101,*),"Number of Theta Elements      : ", NUM_T_ELEMENTS
    WRITE(101,*),"Number of Phi Elements        : ", NUM_P_ELEMENTS
    WRITE(101,*),"-------------------------------------------------"
    WRITE(101,*),"Finite Element Method Degree  : ", DEGREE
    WRITE(101,*),"Spherical Harmonic maximum L  : ", L_LIMIT
    WRITE(101,*),"-------------------------------------------------"
    WRITE(101,*),"Time to Calculate Coefficents : ", InnerTime
    WRITE(101,*),"-------------------------------------------------"


    WRITE(101,*),"L2 Error                      : ", TMP_INTEGRAL
    WRITE(101,*),"Maximum relative error        : ", TMP_MAX
    WRITE(101,*),"Minimum relative error        : ", TMP_MIN

END IF




IF (OUTPUT_PARAMS == 1) THEN

    PRINT*,"For a run with paramters: "
    PRINT*,"-------------------------------------------------"
    PRINT*,"Source Function                 : ",Source_Function_Flag
    PRINT*,"-------------------------------------------------"
    PRINT*,"Number of Radial Elements       : ", NUM_R_ELEMENTS
    PRINT*,"Number of Theta Elements        : ", NUM_T_ELEMENTS
    PRINT*,"Number of Phi Elements          : ", NUM_P_ELEMENTS
    PRINT*,"-------------------------------------------------"
    PRINT*,"Finite Element Method Degree    : ", DEGREE
    PRINT*,"Spherical Harmonic maximum L    : ", L_LIMIT
    PRINT*,"-------------------------------------------------"
    PRINT*,"Time to Calculate Coefficents   : ", InnerTime
    PRINT*,"-------------------------------------------------"


END IF

PRINT*,"L2 Error                        : ", TMP_INTEGRAL
PRINT*,"Maximum relative error          : ", TMP_MAX
PRINT*,"Minimum relative error          : ", TMP_MIN
WRITE(*,'(/ / /)')


END SUBROUTINE L_TWO_ERROR_3D






 !+103+############################################################################!
!                                                                                   !
!   L_INF_ERROR -   Calculates the L_one Error.  Good for analyzing error along     !
!                       a radial line outward.                                      !
!                                                                                   !
 !#################################################################################!
SUBROUTINE  L_INF_ERROR(L_INF_ERR, CUR, ITERS)

INTEGER, INTENT(IN)                                             :: CUR, ITERS
REAL(KIND = idp), DIMENSION(0:ITERS),INTENT(INOUT)              :: L_INF_ERR


INTEGER, PARAMETER                                              :: POWER = 5
INTEGER                                                         :: POINTS, i
REAL(KIND = idp)                                                :: deltar, rloc, theta, phi
REAL(KIND  = idp)                                            :: potential


theta = 0.0_idp
phi = 0.0_idp

POINTS = 2**POWER

deltar = (R_OUTER - R_INNER)/POINTS

L_INF_ERR(CUR) = 0.0_idp


DO i = 1,POINTS

    rloc = (i-1)*deltar+R_INNER

    potential = Calculate_Potential_At_Location(rloc,0.0_idp,0.0_idp)

    L_INF_ERR(CUR) = MAX(L_INF_ERR(CUR), ABS(potential-Analytic_Solution(rloc, theta, phi)))


END DO




END SUBROUTINE L_INF_ERROR








 !+201+############################################################################!
!                                                                                   !
!   MINMAX_RELATIVE_ERROR - Subroutine that finds the maximum and minimum relative  !
!                               error on entire computational domain.               !
!                                                                                   !
 !#################################################################################!
SUBROUTINE MINMAX_RELATIVE_ERROR()


REAL(KIND = idp)                                    ::  TMP_MAX,                        &
                                                        TMP_MIN









END SUBROUTINE MINMAX_RELATIVE_ERROR


END MODULE Error_Analysis_Functions
