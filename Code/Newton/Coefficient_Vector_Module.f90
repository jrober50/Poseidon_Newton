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
            ONLY :  NUM_R_ELEMENTS, NUM_R_NODES,                                    &
                    DEGREE, L_LIMIT, MATRIX_FORMAT, LINEAR_SOLVER,                  &
                    Source_Vector,                                                  &
                    STF_MAT, STF_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, STF_NNZ,       &
                    Coefficient_Vector,                                             &
                    Matrix_Cholesky_Factorized_Flag


USE Boundary_Conditions_Module, &
            ONLY :  DIRICHLET_BC, NEUMANN_BC, DIRICHLET_BC_CCS, NEUMANN_BC_CCS,     &
                    DIRICHLET_BC_CHOL


USE Linear_Solvers_And_Preconditioners, &
            ONLY :  PRECOND_CONJ_GRAD_FULL, PRECOND_CONJ_GRAD_CCS



USE Cholesky_Module,    &
            ONLY :  Cholesky_Factorization, CCS_Forward_Substitution, CCS_Back_Substitution



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

INTEGER                                                                     ::  INFO, LDAB
INTEGER, DIMENSION(NUM_R_NODES)                                             ::  IPIV


REAL(KIND = idp)                                                            ::  SCALE_FACTOR
COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1)                             ::  WORK_VEC, WORK_VECB



REAL(KIND = idp), ALLOCATABLE, DIMENSION(:,:)                               ::  WORK_MAT

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                                 ::  WORK_ELEM_VAL




INTEGER                                                                     ::  NNZ
INTEGER                                                                     ::  l, m, k
INTEGER                                                                     ::  Guess_Flag


REAL(KIND = idp)       :: timea, timeb


IF (MATRIX_FORMAT =='FULL') THEN

    ALLOCATE (WORK_MAT(0:NUM_R_NODES-1, 0:NUM_R_NODES-1))

ELSE IF (MATRIX_FORMAT == 'CCS' ) THEN

    NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1
    ALLOCATE(WORK_ELEM_VAL(0:NNZ-1))

END IF



WORK_VECB = 0









IF (( LINEAR_SOLVER == "CHOL" ) .AND. (Matrix_Cholesky_Factorized_Flag .EQV. .FALSE.) ) THEN


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


        IF (LINEAR_SOLVER =='FULL') THEN
            !####################################!
            !
            !           Full Matrix Solver       !
            !
            !####################################!


        WORK_MAT = STF_MAT(:,:,l)
        WORK_VEC = Source_Vector(:,m,l)





        CALL DIRICHLET_BC(WORK_MAT, WORK_VEC, l, m)

        CALL NEUMANN_BC(l, WORK_VEC)




        CALL DGESV(NUM_R_NODES, 1, WORK_MAT, NUM_R_NODES, IPIV, WORK_VEC, NUM_R_NODES, INFO)
        IF (INFO > 0) THEN
            print*,"DGESV has failed with INFO = ",INFO
        END IF




        !CALL PRECOND_CONJ_GRAD_FULL(WORK_MAT, WORK_VEC)









        ELSE IF (LINEAR_SOLVER == 'CCS') THEN
            !#######################################################################!
            !                                                                       !
            !           CCS Preconditioned Conjugate Gradient Matrix Solver         !
            !                                                                       !
            !#######################################################################!

            WORK_ELEM_VAL = STF_ELEM_VAL(:,l)
            WORK_VEC = Source_Vector(:,m,l)



            CALL DIRICHLET_BC_CCS(  NUM_R_NODES,   NNZ, l, m,                               &
                                    WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, WORK_VEC)



            CALL NEUMANN_BC_CCS(    NUM_R_NODES, NNZ, l, m,                                 &
                                    WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, WORK_VEC)





            GUESS_FLAG = 0  !!! 0 Means no guess, 1 means guess


            IF ( 1 == 0) THEN

                WORK_VECB = WORK_VEC

                CALL PRECOND_CONJ_GRAD_CCS(NUM_R_NODES, NNZ, WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, &
                                            WORK_VECB, GUESS_FLAG, WORK_VECB,0)


                GUESS_FLAG = 1

            END IF




            CALL PRECOND_CONJ_GRAD_CCS(NUM_R_NODES, NNZ, WORK_ELEM_VAL, STF_COL_PTR, STF_ROW_IND,       &
                                        WORK_VEC,GUESS_FLAG, WORK_VECB,0)





        ELSE IF (LINEAR_SOLVER == "CHOL") THEN
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

        END IF



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
