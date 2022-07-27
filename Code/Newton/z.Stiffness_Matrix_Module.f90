   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Stiffness_Matrix_Module                                                      !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions and subroutines associated with the stiffness matrix !##!
!##!        involved in the linear system formed by the expansions. This includes   !##!
!##!        functions that build the matrix in various storage formats.             !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Stiffness_Matrix                                           !##!
!##!    +102+   Deallocate_Stiffness_Matrix                                         !##!
!##!    +103+   Initialize_Stiffness_Matrix_Value                                   !##!
!##!    +104+   Generate_Stiffness_Matrix                                           !##!
!##!                                                                                !##!
!##!    +201+   Generate_Stiffness_Matrix_FULL                                      !##!
!##!    +202+   Generate_Stiffness_Matrix_CCS                                       !##!
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
            ONLY : idp


USE Global_Variables_And_Parameters, &
            ONLY :  NUM_R_ELEMENTS, DEGREE, L_LIMIT, NUM_R_NODES,           &
                    MATRIX_FORMAT,                                          &
                    Stiffness_Matrix_Initialized_Flag,                      &
                    rlocs,                                                  &
                    STF_MAT_Integrals, STF_MAT,                             &   !!! For Full Matrix Storage
                    STF_ELEM_VAL, STF_COL_PTR, STF_ROW_IND, STF_NNZ             !!! For CCS Matrix Storage

USE Additional_Functions_Module, &
            ONLY : Lagrange_Poly, Lagrange_Poly_Deriv, Initialize_LGL_Quadrature




IMPLICIT NONE


CONTAINS
!+101+###########################################################################!
!                                                                                !
!                  ALLOCATE_STF_MAT_MATRIX - Allocate the Source vector          !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Stiffness_Matrix






ALLOCATE(STF_MAT_Integrals(0:3, 0:DEGREE, 0:DEGREE))



!!! Allocate Matrix Space dependent on format !!!
IF (MATRIX_FORMAT =='FULL') THEN





    ALLOCATE(STF_MAT(0:NUM_R_NODES-1, 0:NUM_R_NODES-1, 0:L_LIMIT))






ELSE IF (MATRIX_FORMAT == 'CCS') THEN





    STF_NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1

    !!! Allocate CCS Matrix Arrays !!!
    ALLOCATE(STF_ELEM_VAL(0:STF_NNZ-1,0:L_LIMIT), STF_ROW_IND(0:STF_NNZ-1), STF_COL_PTR(0:NUM_R_NODES))



END IF






END SUBROUTINE Allocate_Stiffness_Matrix






!+102+###########################################################################!
!                                                                                !
!                  DEALLOCATE_STF_MAT_MATRIX - Allocate the Source vector          !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Stiffness_Matrix


DEALLOCATE(STF_MAT_Integrals)



!!! Deallocate Matrix Space dependent on format !!!
IF (MATRIX_FORMAT =='FULL') THEN


    DEALLOCATE(STF_MAT)


ELSE IF (MATRIX_FORMAT == 'CCS') THEN

    !!! Allocate CCS Matrix Arrays !!!
    DEALLOCATE(STF_ELEM_VAL, STF_ROW_IND, STF_COL_PTR)


END IF


END SUBROUTINE Deallocate_Stiffness_Matrix





!+103+##########################################################################!
!                                                                               !
!   Initialize_Stiffness_Matrix - Calculates and returns reusable integrals     !
!       needed to generate the stiffness matrix in any storage format.          !
!                                                                               !
!       *   Run Everytime Num_R_Elements, and/or Degree is set/reset   *        !
!                                                                               !
!###############################################################################!
SUBROUTINE Initialize_Stiffness_Matrix_Values()


INTEGER                                             ::  i, j, k, OrdofQuad, HERE


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP, DLagP, xlocP, weightP
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  xlocQ, weightQ



STF_MAT_Integrals = 0.0_idp


OrdofQuad = DEGREE + 10

ALLOCATE(xlocQ(0:OrdofQuad), weightQ(0:OrdofQuad))

!!! Generate nodes for Lagrange Polynomials !!!
CALL Initialize_LGL_Quadrature(DEGREE, xlocP, weightP)

!!! Generate nodes and weights for Quadrature !!!
CALL Initialize_LGL_Quadrature(OrdofQuad, xlocQ, weightQ)


DO k = 0,OrdofQuad

    !!! Generate Lagrange Polynomials for xlocQ(k) using xlocP for nodes !!!
    LagP = Lagrange_Poly(xlocQ(k), DEGREE, xlocP)

    !!! Generate First Derivatives of Lagrange Polynomials for xlocQ(k) using xlocP for nodes !!!
    DLagP = Lagrange_Poly_Deriv(xlocQ(k), DEGREE, xlocP)

    DO j = 0,DEGREE


        DO i = 0,DEGREE


            !!! Integral of LagP*LagP in x-space. Mapped to r-space in GENERATE_STF_MAT
            STF_MAT_Integrals(0,i,j) = STF_MAT_Integrals(0,i,j) + weightQ(k)*LagP(i)*LagP(j)



            !!! Int of r^2 * DLagP * DLagP dr in x-space.
            !!! Mapping to x-space creates three integral
            !!! 1) Int DLagP*DLagP
            STF_MAT_Integrals(1,i,j) = STF_MAT_Integrals(1,i,j) + weightQ(k)*DLagP(i)*DLagP(j)

            !!! 2) Int x*DLagP*DLagP
            STF_MAT_Integrals(2,i,j) = STF_MAT_Integrals(2,i,j) + weightQ(k)*xlocQ(k)*DLagP(i)*DLagP(j)

            !!! 3) Int x^2*DLagP*DLagP
            STF_MAT_Integrals(3,i,j) = STF_MAT_Integrals(3,i,j) + weightQ(k)*xlocQ(k)*xlocQ(k)*DLagP(i)*DLagP(j)

            !!! In GENERATE_STF_MAT these integrals are mapped back to r-space for given a given range of r.



        END DO

    END DO

END DO





IF (MATRIX_FORMAT == 'CCS') THEN
    !
    !   Because the structure of the stiffness matrix is known, we can
    !       initialize the ROW_IND, and COL_PTR arrays now, and will
    !       fill ELEM_VAL array when GENERATE_STF_MAT_CCS is called.
    !

      !                             !
     !!                             !!
    !!!    COL_PTR INITIALIZATION   !!!
     !!                             !!
      !                             !

    STF_COL_PTR(0) = 0
    STF_COL_PTR(1) = STF_COL_PTR(0) + (DEGREE+1)
    HERE = 2

    DO i = 1,NUM_R_ELEMENTS-1

        DO j = 1,DEGREE - 1

            STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) +  (DEGREE + 1)

            Here = Here + 1

        END DO

        STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) + (2*DEGREE + 1)

        Here = Here + 1

    END DO

    DO i = 1, DEGREE

        STF_COL_PTR(Here) = STF_COL_PTR(Here - 1) + (DEGREE + 1)

        Here = Here + 1

    END DO



      !                             !
     !!                             !!
    !!!    ROW_IND INITIALIZATION   !!!
     !!                             !!
      !                             !

    Here = 0

    DO i = 0, NUM_R_ELEMENTS - 1

        DO j = 0,DEGREE - 1

            DO k = 0, DEGREE

                STF_ROW_IND(Here) = i*DEGREE + k

                Here = Here + 1

            END DO

        END DO

        DO k = 0,DEGREE - 1

            STF_ROW_IND(Here) = i*DEGREE + k

            Here = Here + 1

        END DO

    END DO


    STF_ROW_IND(Here) = DEGREE * NUM_R_ELEMENTS


END IF


Stiffness_Matrix_Initialized_Flag = .TRUE.




END SUBROUTINE Initialize_Stiffness_Matrix_Values








!+104+###########################################################################!
!                                                                                !
!                  Generate_Stiffness_Matrix - Allocate the Source vector          !
!                                                                                !
!################################################################################!
SUBROUTINE Generate_Stiffness_Matrix





!!! Deallocate Matrix Space dependent on format !!!
IF (MATRIX_FORMAT =='FULL') THEN


    CALL Generate_Stiffness_Matrix_Full()


ELSE IF (MATRIX_FORMAT == 'CCS') THEN


    CALL Generate_Stiffness_Matrix_CCS()


END IF


END SUBROUTINE Generate_Stiffness_Matrix















!+201+#################################################################!
!                                                                       !
!     Generate_Stiffness_Matrix_Full -  Generates the STF_MAT Matrix    !
!                                                                       !
!#######################################################################!
SUBROUTINE Generate_Stiffness_Matrix_FULL()


INTEGER                     ::  i, j, iOrd, i1Ord

REAL(KIND = idp)                                ::  deltar, deltarb, rplusr
REAL(KIND = idp), DIMENSION(0:DEGREE)            ::  xlocP, weightP, LagP

STF_MAT = 0.0_idp


DO j = 0,L_LIMIT
    DO i = 0,NUM_R_ELEMENTS-1

        deltar = rlocs(i+1) - rlocs(i)
        rplusr = rlocs(i+1) + rlocs(i)


        iOrd = i*DEGREE
        i1Ord = (i+1)*DEGREE


        STF_MAT(iOrd:i1Ord, iOrd:iOrd, j) = STF_MAT(iOrd:i1Ord, iOrd:i1Ord, j) &
                                        + j*(j+1)*(deltar/2.0_idp)*STF_MAT_Integrals(0,:,:) &
                                        - ((rplusr*rplusr)/(2.0_idp*deltar)*STF_MAT_Integrals(1,:,:)   &
                                            + rplusr*STF_MAT_Integrals(2,:,:) &
                                            + (deltar/2.0_idp)*STF_MAT_Integrals(3,:,:))




        !!! Cartesian !!!
!        STF_MAT(iOrd:i1Ord, iOrd:i1Ord, j) = STF_MAT(iOrd:i1Ord, iOrd:iOrd, j) &
!                                                - (deltar/2.0_idp)*STF_MAT_Integrals(0,:,:)



    END DO
END DO






END SUBROUTINE Generate_Stiffness_Matrix_FULL







!+202+######################################################################!
!                                                                           !
!   Generate_Stiffness_Matrix_CCS -  Generates the stiffness matrix in the  !
!                                      Compress Column Storage (CCS) format !
!                                                                           !
!###########################################################################!
SUBROUTINE Generate_Stiffness_Matrix_CCS()


INTEGER                     ::  i, j, h, k, Here

REAL(KIND = idp)                                ::  deltar, deltarb, rplusr
REAL(KIND = idp), DIMENSION(0:DEGREE)            ::  xlocP, weightP, LagP



STF_ELEM_VAL = 0.0_idp

Here = 0
DO j = 0,L_LIMIT

    Here = 0
    DO i = 0,NUM_R_ELEMENTS-1

        deltar = rlocs(i+1) - rlocs(i)
        rplusr = rlocs(i+1) + rlocs(i)


        DO h = 0, DEGREE


            DO k = 0, DEGREE



                    STF_ELEM_VAL(Here, j)= STF_ELEM_VAL(Here, j) &
                                                + j*(j+1)*(deltar/2.0_idp)*STF_MAT_Integrals(0,h,k) &
                                                + (     (rplusr*rplusr)/(2.0_idp*deltar)            &
                                                            * STF_MAT_Integrals(1,h,k)              &
                                                        + rplusr*STF_MAT_Integrals(2,h,k)           &
                                                        + (deltar/2.0_idp)*STF_MAT_Integrals(3,h,k) )

                    Here = Here + 1



            END DO

        END DO

        Here = Here - 1     !!! Take one step back due to overlap of first value
                            !!! in next element with last in current element

    END DO

END DO





END SUBROUTINE Generate_Stiffness_Matrix_CCS

































END MODULE Stiffness_Matrix_Module