   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE IO_Functions_Module                                                          !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the subroutines used to facilitate the output or input of parts    !##!
!##!        of the linear system to/from a file.                                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   OUTPUT_COEFFS_TO_FILE                                               !##!
!##!    +102+   INPUT_COEFFS_FROM_FILE                                              !##!
!##!    +103+   OUTPUT_SRC_TO_FILE                                                  !##!
!##!    +104+   OUTPUT_STF_MATRIX_TO_FILE                                           !##!
!##!                                                                                !##!
!##!    +201+   MATRIX_TO_AIJ                                                       !##!
!##!    +202+   MATRIX_TO_CRS                                                       !##!
!##!    +203+   MATRIX_TO_CCS                                                       !##!
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
                    ONLY :  NUM_R_NODES, NUM_R_ELEMENTS, DEGREE, L_LIMIT,   &
                            Source_Vector,                                  &
                            STF_MAT, STF_ELEM_VAL,                          &
                            Coefficient_Vector

USE COMMON_IO, &
                    ONLY : FILE_DATA_TYPE

IMPLICIT NONE


!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS







!+101+######################################################################!
!                                                                           !
!   OUTPUT_COEFFS_TO_FILE - Write coefficents to file.  Useful for testing  !
!                                                                           !
!###########################################################################!
SUBROUTINE OUTPUT_COEFFS_TO_FILE(FILE)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE



INTEGER                                                     ::  unit
INTEGER                                                     ::  i, j, k


unit = FILE%FILE_UNIT

!WRITE(unit,*) NUM_R_NODES, L_LIMIT


DO i = 0,L_LIMIT

    DO j = -L_LIMIT,L_LIMIT


        WRITE(unit,*) i, j, REALPART(Coefficient_Vector(:,j,i))



    END DO ! j Loop

END DO ! i Loop





END SUBROUTINE OUTPUT_COEFFS_TO_FILE











!+102+######################################################################!
!                                                                           !
!   INPUT_COEFFS_FROM_FILE - Read coefficents from file. Useful for testing !
!                                                                           !
!###########################################################################!
SUBROUTINE INPUT_COEFFS_FROM_FILE(FILE, COEFFS_TMP)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE
COMPLEX(KIND = idp), INTENT(INOUT), DIMENSION(0:NUM_R_NODES-1, -L_LIMIT:L_LIMIT, 0:L_LIMIT) :: COEFFS_TMP



INTEGER                                                     ::  unit
INTEGER                                                     ::  i, j, k

INTEGER                                                     :: NUM_R_NODES_TMP, L_LIMIT_TMP





unit = FILE%FILE_UNIT

READ(unit,*) NUM_R_NODES_TMP, L_LIMIT_TMP


DO i = 0,L_LIMIT

    DO j = -L_LIMIT,L_LIMIT


        READ(unit,*) COEFFS_TMP(:,j,i)



    END DO ! j Loop

END DO ! i Loop





END SUBROUTINE INPUT_COEFFS_FROM_FILE














!+103+#######################################################################!
!                                                                           !
!   OUTPUT_SRC_TO_FILE - Write coefficents to file.  Useful for testing  !
!                                                                           !
!###########################################################################!
SUBROUTINE OUTPUT_SRC_TO_FILE(FILE)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE



INTEGER                                                     ::  unit
INTEGER                                                     ::  i, j, k


unit = FILE%FILE_UNIT

!WRITE(unit,*) NUM_R_NODES, L_LIMIT


DO i = 0,L_LIMIT

    DO j = -L_LIMIT,L_LIMIT

        DO k = 0,NUM_R_NODES-1


            WRITE(unit,*),i,j,k, REALPART(Source_Vector(k,j,i)), IMAGPART(Source_Vector(k,j,i))

        END DO ! k Loop

    END DO ! j Loop

END DO ! i Loop





END SUBROUTINE OUTPUT_SRC_TO_FILE





!+104+#######################################################################!
!                                                                           !
!   OUTPUT_STF_MATRIX_TO_FILE - Write coefficents to file.  Useful for testing  !
!                                                                           !
!###########################################################################!
SUBROUTINE OUTPUT_STF_MATRIX_TO_FILE(FILE)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE



INTEGER                                                     ::  unit
INTEGER                                                     ::  i, j, k, NNZ


unit = FILE%FILE_UNIT

!WRITE(unit,*) NUM_R_NODES, L_LIMIT

NNZ = NUM_R_ELEMENTS*(DEGREE + 1)*(DEGREE + 1) - NUM_R_ELEMENTS + 1



DO i = 0,L_LIMIT

    DO j = 0,NNZ-1


        WRITE(unit,*) STF_ELEM_VAL(j,i)



    END DO ! j Loop

END DO ! i Loop





END SUBROUTINE OUTPUT_STF_MATRIX_TO_FILE









!+201+#################################################################################!

!                  MATRIX_TO_AIJ - Print matrix to file in AIJ format                   !

!#######################################################################################!
SUBROUTINE MATRIX_TO_AIJ(WORK_MAT, M, N, NAME)


INTEGER, INTENT(IN)                                         :: M, N

REAL(KIND = idp), DIMENSION(0:M-1,0:N-1), INTENT(IN)        :: WORK_MAT


CHARACTER(LEN = 30),INTENT(IN)                              :: NAME


INTEGER                                                     :: i, j, here
INTEGER                                                     :: NUM_NNZ, NUM_ROW_NNZ, NUM_COL_NNZ
LOGICAL, DIMENSION(0:M-1,0:N-1)                             :: NNZ_MAT

INTEGER, ALLOCATABLE, DIMENSION(:,:)                        :: LOC
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 :: ELEM_VAL


!!! Create Map of non-zeros using Logical Array !!!
NNZ_MAT = (WORK_MAT .NE. 0)

!!! Count non-zeros using non-zero map !!!
NUM_NNZ = COUNT(NNZ_MAT)


ALLOCATE(LOC(0:1,0:NUM_NNZ-1), ELEM_VAL(0:NUM_NNZ-1))

here = 0
DO i = 0,M-1
    DO j = 0,N-1

        IF (NNZ_MAT(j,i) .EQV. .TRUE.) THEN

            LOC(0,here) = i
            LOC(1,here) = j
            ELEM_VAL(here) = WORK_MAT(j,i)

            here = here + 1

        END IF

    END DO
END DO



OPEN(unit = 42, file = NAME)
WRITE(42,*), M, N, NUM_NNZ

DO i = 0, NUM_NNZ-1

    WRITE(42,*),LOC(0,i),LOC(1,i),ELEM_VAL(i)

END DO

CLOSE(unit = 42)
PRINT*,"Matrix in AIJ format is at ",NAME



END SUBROUTINE MATRIX_TO_AIJ









!+202+###################################################################################!

!                  MATRIX_TO_CRS - Print matrix in CRS format to two files              !

!                      First File: COL_ELEM.CSR - contains column indicator and element !
!                                                    value arrays                       !

!                      Second File:  ROW_PTR.CSR - contains row pointer array           !

!#######################################################################################!
SUBROUTINE MATRIX_TO_CRS(WORK_MAT, M, N)


INTEGER, INTENT(IN)                                         :: M, N

REAL(KIND = idp), DIMENSION(0:M-1,0:N-1), INTENT(IN)        :: WORK_MAT


CHARACTER(LEN = 20)                                         :: filenameA = 'OUTPUT/COL_ELEM.CRS'
CHARACTER(LEN = 20)                                         :: filenameB = 'OUTPUT/ROW_PTR.CRS'


INTEGER                                                     :: i, j, here
INTEGER                                                     :: NUM_NNZ, NUM_ROW_NNZ, NUM_COL_NNZ
LOGICAL, DIMENSION(0:M-1,0:N-1)                             :: NNZ_MAT

INTEGER, ALLOCATABLE, DIMENSION(:)                          :: ROW_PTR, COL_IND
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 :: ELEM_VAL



!!! Use Logical array to create non-zero map !!!
NNZ_MAT = (WORK_MAT .NE. 0)


!!! Use non-zero map to count non-zeros !!!
NUM_NNZ = COUNT(NNZ_MAT)



!NUM_ROW_NNZ = COUNT(NNZ_MAT(0,:))
!NUM_COL_NNZ = COUNT(NNZ_MAT(:,1))

ALLOCATE(COL_IND(0:NUM_NNZ-1), ELEM_VAL(0:NUM_NNZ-1), ROW_PTR(0:M))

here = 0
DO i = 0,M-1

    ROW_PTR(i) = here

    NUM_ROW_NNZ = COUNT(NNZ_MAT(i,:))

    IF (NUM_ROW_NNZ .NE. 0) THEN

        DO j = 0,M

            IF (NNZ_MAT(j,i) .EQV. .TRUE.) THEN

                COL_IND(here) = j
                ELEM_VAL(here) = WORK_MAT(j,i)

                here = here + 1

            END IF

        END DO


    END IF


END DO
ROW_PTR(M) = NUM_NNZ



!!! Output Column indicator and Element Values file !!!
OPEN(unit = 42, file = filenameA)

DO i = 0, NUM_NNZ-1

    WRITE(42,*),COL_IND(i),ELEM_VAL(i)

END DO

CLOSE(unit = 42)
PRINT*,"Matrix in CRS format is at ",filenameA





!!! Output row pointer file  !!!
OPEN(unit = 42, file = filenameB)
DO i = 0,M

    WRITE(42,*), ROW_PTR(i)

END DO
CLOSE(unit = 42)
print*," and ", filenameB





END SUBROUTINE MATRIX_TO_CRS
















!+203+##################################################################################!

!                  MATRIX_TO_CCS - Print matrix in CCS format to two files              !

!                      First File: ROW_ELEM.CSR - contains row indicator and element    !
!                                                    value arrays                       !

!                      Second File:  COL_PTR.CSR - contains column pointer array        !

!#######################################################################################!
SUBROUTINE MATRIX_TO_CCS(WORK_MAT, M, N)


INTEGER, INTENT(IN)                                         :: M, N

REAL(KIND = idp), DIMENSION(0:M-1,0:N-1), INTENT(IN)        :: WORK_MAT


CHARACTER(LEN = 20)                                         :: filenameA = 'OUTPUT/ROW_ELEM.CCS'
CHARACTER(LEN = 20)                                         :: filenameB = 'OUTPUT/COL_PTR.CCS'


INTEGER                                                     :: i, j, here
INTEGER                                                     :: NUM_NNZ, NUM_ROW_NNZ, NUM_COL_NNZ
LOGICAL, DIMENSION(0:M-1,0:N-1)                             :: NNZ_MAT

INTEGER, ALLOCATABLE, DIMENSION(:)                          :: COL_PTR, ROW_IND
REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)                 :: ELEM_VAL



!!! Use Logical array to create non-zero map !!!
NNZ_MAT = (WORK_MAT .NE. 0)


!!! Use non-zero map to count non-zeros !!!
NUM_NNZ = COUNT(NNZ_MAT)



!NUM_ROW_NNZ = COUNT(NNZ_MAT(0,:))
!NUM_COL_NNZ = COUNT(NNZ_MAT(:,1))

ALLOCATE(ROW_IND(0:NUM_NNZ-1), ELEM_VAL(0:NUM_NNZ-1), COL_PTR(0:M))

here = 0
DO i = 0,M-1

    COL_PTR(i) = here

    NUM_COL_NNZ = COUNT(NNZ_MAT(:,i))

    IF (NUM_COL_NNZ .NE. 0) THEN

        DO j = 0,M

            IF (NNZ_MAT(j,i) .EQV. .TRUE.) THEN

                ROW_IND(here) = j
                ELEM_VAL(here) = WORK_MAT(j,i)

                here = here + 1

            END IF

        END DO




    END IF


END DO
COL_PTR(N) = NUM_NNZ



!!! Output Column indicator and Element Values file !!!
OPEN(unit = 42, file = filenameA)

DO i = 0, NUM_NNZ-1

    WRITE(42,*),ROW_IND(i),ELEM_VAL(i)

END DO

CLOSE(unit = 42)
PRINT*,"Matrix in CCS format is at ",filenameA






!!! Output row pointer file  !!!
OPEN(unit = 42, file = filenameB)
DO i = 0,M

    WRITE(42,*), COL_PTR(i)

END DO
CLOSE(unit = 42)
print*," and ", filenameB




END SUBROUTINE MATRIX_TO_CCS














!########################################################################################!
!
!           SCCS_TO_FULL - take matrix in CCS format and creates full matrix
!
!########################################################################################!
SUBROUTINE SCCS_TO_FULL(NUM_COL, NUM_ROW, NUM_NZELEM, ROW_IND, COL_PTR, ELEM_VAL, A)


INTEGER, INTENT(IN)                 :: NUM_COL, NUM_ROW, NUM_NZELEM


INTEGER, DIMENSION(0:NUM_COL), INTENT(IN)                           :: COL_PTR
INTEGER, DIMENSION(0:NUM_NZELEM-1), INTENT(IN)                      :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NUM_NZELEM-1), INTENT(IN)             :: ELEM_VAL

REAL(KIND = idp), DIMENSION(0:NUM_ROW-1,0:NUM_COL-1), INTENT(INOUT) :: A


INTEGER                                                         :: i, j




A = 0.0_idp

DO i = 0, NUM_COL-1

    DO j = COL_PTR(i),COL_PTR(i+1)-1

        A(ROW_IND(j),i) = ELEM_VAL(j)
        A(i,ROW_IND(j)) = ELEM_VAL(j)

    END DO

END DO




END SUBROUTINE SCCS_TO_FULL







!########################################################################################!
!
!           CCS_TO_FULL - take matrix in CCS format and creates full matrix
!
!########################################################################################!
SUBROUTINE CCS_TO_FULL(NUM_COL, NUM_ROW, NUM_NZELEM, ROW_IND, COL_PTR, ELEM_VAL, A)


INTEGER, INTENT(IN)                 :: NUM_COL, NUM_ROW, NUM_NZELEM


INTEGER, DIMENSION(0:NUM_COL), INTENT(IN)                           :: COL_PTR
INTEGER, DIMENSION(0:NUM_NZELEM-1), INTENT(IN)                      :: ROW_IND
REAL(KIND = idp), DIMENSION(0:NUM_NZELEM-1), INTENT(IN)             :: ELEM_VAL

REAL(KIND = idp), DIMENSION(0:NUM_ROW-1,0:NUM_COL-1), INTENT(INOUT) :: A


INTEGER                                                             :: i, j



A = 0.0_idp

DO i = 0, NUM_COL-1

    DO j = COL_PTR(i),COL_PTR(i+1)-1

        !PRINT*,ROW_IND(j), i, ELEM_VAL(j)
        A(ROW_IND(j),i) = ELEM_VAL(j)

    END DO

END DO




END SUBROUTINE CCS_TO_FULL






END MODULE IO_Functions_Module