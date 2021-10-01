   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Coefficient_Vector_Module                                                    !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions and subroutines that concern the coefficient vector, !##!
!##!        most importantly, the main routine responsible for solving the linear   !##!
!##!        system of equations.                                                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Coefficient_Vector                                         !##!
!##!    +102+   Deallocate_Coefficeint_Vector                                       !##!
!##!                                                                                !##!
!##!    +201+   Calculate_Coefficient_Vector                                        !##!
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
            ONLY :  idp



USE Global_Variables_And_Parameters, &
            ONLY :  NUM_R_ELEMENTS, NUM_R_NODES,                            &
                    DEGREE, L_LIMIT, MATRIX_FORMAT, LINEAR_SOLVER,          &
                    Source_Vector,                                          &
                    STF_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, STF_NNZ,        &
                    Coefficient_Vector,                                     &
                    Matrix_Cholesky_Factorized_Flag


USE Boundary_Conditions_Module, &
            ONLY :  NEUMANN_BC_CCS,             &
                    DIRICHLET_BC_CHOL





USE Cholesky_Module,    &
            ONLY :  Cholesky_Factorization

USE CCS_Forward_Substitution_Module, &
            ONLY :  CCS_Forward_Substitution

USE CCS_Back_Substitution_Module, &
            ONLY :  CCS_Back_Substitution



IMPLICIT NONE



CONTAINS
!+101+############################################################################!
!                                                                                !
!               ALLOCATE_COEFFS - Allocate the coefficent vectors                !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Coefficient_Vector()

ALLOCATE( Coefficient_Vector(0:NUM_R_NODES-1, -L_LIMIT:L_LIMIT, 0:L_LIMIT))


END SUBROUTINE Allocate_Coefficient_Vector



!+102+###########################################################################!
!                                                                                !
!             DEALLOCATE_COEFFS - Allocate the coefficent vectors                !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Coefficient_Vector()


DEALLOCATE( Coefficient_Vector )


END SUBROUTINE Deallocate_Coefficient_Vector









!+201+##########################################################################!
!                                                                               !
!      Calculate_Coefficient_Vector - Solve linear system STF_MAT*Coeffs = Src  !
!                                                                               !
!###############################################################################!
SUBROUTINE Calculate_Coefficient_Vector()

INTEGER                                                         ::  INFO, LDAB
INTEGER, DIMENSION(NUM_R_NODES)                                 ::  IPIV


REAL(KIND = idp)                                                ::  SCALE_FACTOR
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                 ::  WORK_VEC



REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                   ::  WORK_MAT
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                     ::  WORK_ELEM_VAL



INTEGER                                                         ::  NNZ
INTEGER                                                         ::  l, m, k
INTEGER                                                         ::  Guess_Flag


REAL(KIND = idp)       :: timea, timeb



NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))





IF ( Matrix_Cholesky_Factorized_Flag .EQV. .FALSE. ) THEN


    !
    !   This only needs to be done everytime the radial mesh is defined/redefined.
    !   This performs Cholesky factorization on the stiffness matrix and overwrites
    !   the stiffness matrix variables STF_ELEM_VAL, STF_ROW_IND, and STF_COL_PTR to
    !   represent the factorization matrix, L.  This matrix can then be reused to 
    !   solve the linear system using forward and backward substitution.
    !

    CALL Cholesky_Factorization()
    Matrix_Cholesky_Factorized_Flag = .TRUE.


END IF











DO l = 0,L_LIMIT

    DO m = -l,l


       

        !#######################################################################!
        !                                                                       !
        !               CCS Cholesky Factorization Matrix Solver                !
        !                                                                       !
        !#######################################################################!



        WORK_VEC = Source_Vector(:,m,l)


        CALL DIRICHLET_BC_CHOL( NUM_R_NODES,    &
                                STF_NNZ,        &
                                l, m,           &
                                STF_COL_PTR,    &
                                STF_ROW_IND,    &
                                WORK_VEC        )


        CALL NEUMANN_BC_CCS(    NUM_R_NODES,    &
                                STF_NNZ,        &
                                l, m,           &
                                WORK_ELEM_VAL,  &
                                STF_COL_PTR,    &
                                STF_ROW_IND,    &
                                WORK_VEC        )



        CALL CCS_Forward_Substitution(  NUM_R_NODES,        &
                                        STF_NNZ,            &
                                        STF_ELEM_VAL(:,l),  &
                                        STF_COL_PTR,        &
                                        STF_ROW_IND,        &
                                        WORK_VEC            )


        CALL CCS_Back_Substitution( NUM_R_NODES,        &
                                    STF_NNZ,            &
                                    STF_ELEM_VAL(:,l),  &
                                    STF_COL_PTR,        &
                                    STF_ROW_IND,        &
                                    WORK_VEC            )




        Do k = 0,NUM_R_NODES - 1
            Coefficient_Vector(k,m,l) = WORK_VEC(k)
        END DO



    END DO
END DO





IF (MATRIX_FORMAT =='FULL') THEN

    DEALLOCATE (WORK_MAT)

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    DEALLOCATE(WORK_ELEM_VAL)

END IF




END SUBROUTINE Calculate_Coefficient_Vector












END MODULE Coefficient_Vector_Module
