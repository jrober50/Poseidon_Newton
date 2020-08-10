   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Boundary_Conditions_Module                                                   !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!   Contains the subroutines used to calculate and impliment boundary conditions !##!
!##!    on the system being solved.                                                 !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   DIRICHLET_BC                                                        !##!
!##!    +101+   NEUMANN_BC                                                          !##!
!##!                                                                                !##!
!##!    +201+   DIRICHLET_BC_CSS                                                    !##!
!##!    +202+   NEUMANN_BC_CSS                                                      !##!
!##!                                                                                !##!
!##!    +301+   DIRICHLET_BC_CHOL                                                   !##!
!##!                                                                                !##!
!##!    +401+   BC_Integral                                                         !##!
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
                ONLY :  NUM_R_NODES, R_INNER, R_OUTER, DEGREE,                          &
                        INNER_BC_TYPE, OUTER_BC_TYPE,                                   &
                        INNER_BC_SET_FLAG, OUTER_BC_SET_FLAG,                           &
                        INNER_DIR_BC_INPUT, INNER_NEU_BC_INPUT,                         &
                        OUTER_DIR_BC_INPUT, OUTER_NEU_BC_INPUT,                         &
                        INNER_UNIFORM_DIR_BC_FLAG, OUTER_UNIFORM_DIR_BC_FLAG,           &
                        L_LIMIT, NUM_R_ELEMENTS, POWER_A, RHO_O, Source_Function_Flag,  &
                        First_Column_Storage, Last_Column_Storage,                      &
                        Analytic_Solution,                                              &
                        STF_MAT,                                                        &
                        Test_Run_Flag


USE Additional_Functions_Module, &
                ONLY :  Legendre_Poly, Lagrange_Poly, Spherical_Harmonic,               &
                        Norm_Factor, POWER, Map_From_X_Space,                           &
                        Initialize_LG_Quadrature, Initialize_LGL_Quadrature


IMPLICIT NONE

!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS




 !+101+############################################################################!
!                                                                                   !
!       Dirichlet_BC                                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!

SUBROUTINE DIRICHLET_BC(WORK_MAT, WORK_VEC, L, M)

INTEGER, INTENT(IN)                                                                     :: L, M

!REAL(KIND = idp), DIMENSION(0:NUM_R_NODES - 1), INTENT(INOUT)                           :: WORK_VEC
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES - 1), INTENT(INOUT)                           :: WORK_VEC

REAL(KIND = idp), DIMENSION(0:NUM_R_NODES - 1,0:NUM_R_NODES - 1), INTENT(INOUT)         :: WORK_MAT




INTEGER                 :: i, shift


COMPLEX(KIND = idp)                                                         :: BC_Value


IF (INNER_BC_TYPE == "DRCH") THEN


    IF (    Test_Run_Flag .EQV. .TRUE.  )   THEN


        BC_Value = BC_Integral(R_INNER, L, M)


    ELSE

        BC_Value =  INNER_DIR_BC_INPUT

    END IF




    WORK_VEC(0) = BC_Value
    DO i = 1,DEGREE


        WORK_VEC(i) = WORK_VEC(i) - STF_MAT(0,i,L)*BC_Value


    END DO




    WORK_MAT(0,:) = 0.0_idp
    WORK_MAT(:,0) = 0.0_idp
    WORK_MAT(0,0) = 1.0_idp

END IF




shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (OUTER_BC_TYPE  == "DRCH") THEN



    IF (    Test_Run_Flag .EQV. .TRUE.  )   THEN


        BC_Value = BC_Integral(R_INNER, L, M)


    ELSE

        BC_Value =  OUTER_DIR_BC_INPUT

    END IF




    WORK_VEC(NUM_R_NODES - 1) = BC_Value

    DO i = 1,DEGREE-shift


        WORK_VEC(NUM_R_NODES-1-i) = WORK_VEC(NUM_R_NODES-1-i) - STF_MAT(NUM_R_NODES-1, NUM_R_NODES-1-i,L)*BC_Value

    END DO


    WORK_MAT(NUM_R_NODES-1,:) = 0.0_idp
    WORK_MAT(:,NUM_R_NODES-1) = 0.0_idp
    WORK_MAT(NUM_R_NODES-1,NUM_R_NODES-1) = 1.0_idp




END IF






END SUBROUTINE DIRICHLET_BC




 !+102+############################################################################!
!                                                                                   !
!       Neumann_BC                                                                !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE NEUMANN_BC(L_VALUE, WORK_VEC)

INTEGER,                                            INTENT(IN)                      ::  L_VALUE
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES - 1),  INTENT(INOUT)                   ::  WORK_VEC


INTEGER                 :: i,j




IF (    L_VALUE == 0    )   THEN


    IF (INNER_BC_TYPE == "NEUM") THEN

        WORK_VEC(0) = WORK_VEC(0) + R_INNER*R_INNER*INNER_NEU_BC_INPUT

    END IF


    IF (OUTER_BC_TYPE == "NEUM") THEN

        WORK_VEC(NUM_R_NODES-1) = WORK_VEC(NUM_R_NODES-1) + R_OUTER*R_OUTER*OUTER_NEU_BC_INPUT

    END IF



END IF


END SUBROUTINE NEUMANN_BC















 !+201+############################################################################!
!                                                                                   !
!       Dirichlet_BC_CCS                                                            !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_CCS(N, NNZ, L, M, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC


REAL(KIND = idp), DIMENSION(0:NNZ-1),           INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                                     :: i, shift
COMPLEX(KIND = idp)                                                         :: BC_Value






IF (INNER_BC_TYPE == "DRCH") THEN



    IF (    INNER_UNIFORM_DIR_BC_FLAG .EQV. .TRUE.  )   THEN




          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
        IF (    L == 0  )   THEN


            BC_Value = sqrt(4.0_idp*pi)*INNER_DIR_BC_INPUT

        ELSE

            BC_Value = 0.0_idp

        END IF

    ELSE    ! For a Non-Uniform Boundary Value the integrals will have to be evaluated !
            !   Currently not working.

        BC_Value =  BC_Integral(R_INNER, L, M)

    END IF





    !!! MODIFY SOURCE VECTOR !!!


    WORK_VEC(0) = BC_Value

    DO i = COL_PTR(0)+1,COL_PTR(1)-1

        WORK_VEC(i) = WORK_VEC(i) - ELEM_VAL(i)*BC_Value

    END DO



    !!! MODIFY MATRIX !!!

    ELEM_VAL(0) = 1.0_idp

    DO i = 1,DEGREE

        ELEM_VAL(i) = 0.0_idp
        ELEM_VAL(COL_PTR(i)) = 0.0_idp


    END DO




END IF









shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (OUTER_BC_TYPE  == "DRCH") THEN


    IF (    OUTER_UNIFORM_DIR_BC_FLAG .EQV. .TRUE.  )   THEN



          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
        IF (    L == 0  )   THEN

            BC_Value = sqrt(4.0_idp*pi)*OUTER_DIR_BC_INPUT

        ELSE

            BC_Value = 0.0_idp

        END IF

    ELSE    ! NON-Uniform Boundary Value

        BC_Value =  BC_Integral(R_OUTER, L, M)

    END IF



    !!! MODIFY SRC VECTOR !!!
    WORK_VEC(NUM_R_NODES - 1) = BC_Value

    DO i = 1,DEGREE-shift


        WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                        - ELEM_VAL(COL_PTR(N) - 1 -i)*BC_Value



    END DO



    !!! MODIFY MATRIX !!!
    DO i = 1,DEGREE

        ELEM_VAL(NNZ - 1 - i) = 0.0_idp
        ELEM_VAL(COL_PTR(NUM_R_NODES - i )-1) = 0.0_idp


    END DO

    ELEM_VAL(NNZ-1) = 1.0_idp


END IF










END SUBROUTINE DIRICHLET_BC_CCS


















 !+202+############################################################################!
!                                                                                   !
!       Neumann_BC_CCS                                                              !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE NEUMANN_BC_CCS(N, NNZ, L, M, ELEM_VAL, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC


REAL(KIND = idp), DIMENSION(0:NNZ-1),           INTENT(INOUT)               ::  ELEM_VAL


INTEGER                                                                     ::  i, shift

REAL(KIND = idp)                                                            ::  BC_Enc_Mass,    &
                                                                                Shift_Value







IF (    INNER_BC_TYPE .EQ. "NEUM"   ) THEN



    BC_Enc_Mass = INNER_NEU_BC_INPUT



     !                                                   !
    !!   We are asssuming a uniform value on the shell   !!
    !!   so only the L,M = 0 vector is effected          !!
     !                                                   !
    IF (    L .EQ. 0    ) THEN

        !                           !
        !   Calulcate Shift Value   !
        !                           !

        Shift_Value = sqrt(4.0_idp*pi)*BC_Enc_Mass


        WORK_VEC(0) = WORK_VEC(0) - Shift_Value





    END IF

END IF










IF (    OUTER_BC_TYPE .EQ. "NEUM"   ) THEN


    BC_Enc_Mass = OUTER_NEU_BC_INPUT


     !                                                   !
    !!   We are asssuming a uniform value on the shell   !!
    !!   so only the L,M = 0 vector is effected.         !!
     !                                                   !
    IF (    L .EQ. 0    ) THEN

        !                           !
        !   Calulcate Shift Value   !
        !                           !

        Shift_Value = sqrt(4.0_idp*pi)*BC_Enc_Mass


        WORK_VEC(N-1) = WORK_VEC(N-1) + Shift_Value





    END IF

END IF









END SUBROUTINE NEUMANN_BC_CCS













 !+301+############################################################################!
!                                                                                   !
!       Dirichlet_BC_CHOL                                                           !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
 !#################################################################################!
SUBROUTINE DIRICHLET_BC_CHOL(N, NNZ, L, M, COL_PTR, ROW_IND, WORK_VEC)


INTEGER,                                        INTENT(IN)                  ::  N, NNZ, L, M

INTEGER, DIMENSION(0:N),                        INTENT(IN)                  ::  COL_PTR
INTEGER, DIMENSION(0:NNZ-1),                    INTENT(IN)                  ::  ROW_IND

COMPLEX(KIND = idp), DIMENSION(0:N - 1),        INTENT(INOUT)               ::  WORK_VEC




INTEGER                                                                     :: i, shift
COMPLEX(KIND = idp)                                                         :: BC_Value






IF (INNER_BC_TYPE == "DRCH") THEN



    IF (    INNER_UNIFORM_DIR_BC_FLAG .EQV. .TRUE.  )   THEN




          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
        IF (    L == 0  )   THEN


            BC_Value = sqrt(4.0_idp*pi)*INNER_DIR_BC_INPUT

        ELSE

            BC_Value = 0.0_idp

        END IF

    ELSE    ! For a Non-Uniform Boundary Value the integrals will have to be evaluated !
            !   Currently not working.

        BC_Value =  BC_Integral(R_INNER, L, M)

    END IF





    !!! MODIFY SOURCE VECTOR !!!


    WORK_VEC(0) = BC_Value

    DO i = 1,DEGREE

        WORK_VEC(i) = WORK_VEC(i) - First_Column_Storage(i , L)*BC_Value

    END DO



    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




END IF









shift = 0
IF (NUM_R_ELEMENTS .EQ. 1 ) THEN
    shift = 1
END IF

IF (OUTER_BC_TYPE  == "DRCH") THEN


    IF (    OUTER_UNIFORM_DIR_BC_FLAG .EQV. .TRUE.  )   THEN



          !                                                                                       !
         !!   For a uniform boundary condition on a sphere the only spherical harmonic expansion  !!
        !!!   coefficient to be effected will be the l,m = 0,0 coefficient due to the symmetry    !!!
        !!!   of the condition, and the orthogonality of the spherical harmonics functions.       !!!
         !!   The integral over the theta, and phi also produces a sqrt(4*pi) scalling factor.    !!
          !                                                                                       !
        IF (    L == 0  )   THEN

            BC_Value = sqrt(4.0_idp*pi)*OUTER_DIR_BC_INPUT

        ELSE

            BC_Value = 0.0_idp

        END IF

    ELSE    ! NON-Uniform Boundary Value

        BC_Value =  BC_Integral(R_OUTER, L, M)

    END IF



    !!! MODIFY SRC VECTOR !!!
    WORK_VEC(NUM_R_NODES - 1) = BC_Value


    DO i = DEGREE-shift,1,-1


        WORK_VEC(NUM_R_NODES - 1 - i) = WORK_VEC(NUM_R_NODES - 1 - i)                           &
                                        - Last_Column_Storage(i,L)*BC_Value



    END DO





    !!! MODIFY MATRIX !!!

    !!! ALREADY DONE IN CHOLESKY FACTORIZATION !!!




END IF










END SUBROUTINE DIRICHLET_BC_CHOL





















 !+401+####################################################################!
!                                                                           !
!       BC_Integral - Boundary Condition Integral
!                                                                           !
 !#########################################################################!
FUNCTION BC_Integral(R, L, M)


COMPLEX(KIND = idp)                                :: BC_Integral



REAL(KIND = idp), INTENT(IN)                    :: R
INTEGER, INTENT(IN)                             :: L, M

INTEGER                                         :: te, td, pe, pd


INTEGER                                         ::  T_Degree, P_Degree

INTEGER                                         ::  T_SUB_ELEMENTS, P_SUB_ELEMENTS


REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)     ::  P_xlocs, P_locs, P_weights, &
                                                    T_xlocs, T_locs, T_weights, &
                                                    Sub_tlocs, Sub_plocs

REAL(KIND = idp)                                ::  deltasubt, deltasubp, deltat, deltap


COMPLEX(KIND = idp)                             ::  TMP_INTEGRAL_VALUE


T_SUB_ELEMENTS = 32
P_SUB_ELEMENTS = 32

T_Degree = 5
P_Degree = 5





ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))
ALLOCATE(Sub_tlocs(0:T_SUB_ELEMENTS),Sub_plocs(0:P_SUB_ELEMENTS))

CALL Initialize_LGL_Quadrature(P_Degree, P_xlocs, P_weights)
CALL Initialize_LGL_Quadrature(T_Degree, T_xlocs, T_weights)

deltasubt = (pi)/REAL(T_SUB_ELEMENTS)
deltasubp = (2*pi)/REAL(P_SUB_ELEMENTS)

DO te = 0, T_SUB_ELEMENTS
    Sub_tlocs(te) = 0.0_idp + te*deltasubt
END DO



DO pe = 0,P_SUB_ELEMENTS
    Sub_plocs(pe) = 0.0_idp + pe*deltasubp
END DO




TMP_INTEGRAL_VALUE = 0.0_idp
DO te = 0, T_SUB_ELEMENTS - 1

    T_locs = Map_From_X_Space(Sub_tlocs(te), Sub_tlocs(te+1), T_xlocs)
    deltat = Sub_tlocs(te+1)-Sub_tlocs(te)

    DO pe = 0, P_SUB_ELEMENTS - 1

        P_locs = Map_From_X_Space(Sub_plocs(pe), Sub_plocs(pe+1), P_xlocs)
        deltap = Sub_plocs(pe+1) - Sub_plocs(pe)

        DO td = 0, T_Degree

            DO pd = 0, P_Degree

                TMP_INTEGRAL_VALUE = TMP_INTEGRAL_VALUE                                                     &
                                    + T_weights(td)*P_weights(pd)                                           &
                                    * POWER(-1.0_idp,M)*Spherical_Harmonic(L, -M, T_locs(td), P_locs(pd))   &
                                    * sin(T_locs(td))*Analytic_Solution(R, T_locs(td), P_locs(pd))          &
                                    * deltat/2.0_idp*deltap/2.0_idp





            END DO

        END DO

    END DO

END DO




BC_Integral = TMP_INTEGRAL_VALUE




END FUNCTION BC_Integral






END MODULE Boundary_Conditions_Module