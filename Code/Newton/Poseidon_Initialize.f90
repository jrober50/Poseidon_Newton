   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE  Poseidon_Initialize_Module                                           !##!
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



USE Global_Variables_And_Parameters, &
                                ONLY :  DEGREE, L_LIMIT, LM_Length,                     &
                                        R_INNER, R_OUTER,                               &
                                        NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                                        NUM_R_NODES,                                    &
                                        rlocs, tlocs, plocs,                            &
                                        Coefficient_Vector,                             &
                                        Source_Term_Coefficients,                       &
                                        Source_Terms,                                   &
                                        Source_Degrees,                                 &
                                        Num_Quad_DOF,                                   &
                                        Test_Source_Input,                              &
                                        Test_Space_Allocated_Flag,                      &
                                        Stiffness_Matrix_Initialized_Flag,              &
                                        INNER_BC_SET_FLAG, OUTER_BC_SET_FLAG,           &
                                        INNER_BC_TYPE, OUTER_BC_TYPE,                   &
                                        INNER_DIR_BC_INPUT, INNER_NEU_BC_INPUT,         &
                                        OUTER_DIR_BC_INPUT, OUTER_NEU_BC_INPUT,         &
                                        INNER_UNIFORM_DIR_BC_FLAG, OUTER_UNIFORM_DIR_BC_FLAG,   &
                                        RADIAL_MESH_SET_FLAG,                           &
                                        THETA_MESH_SET_FLAG,                            &
                                        PHI_MESH_SET_FLAG,                              &
                                        TEST_RUN_FLAG,                                  &
                                        First_Column_Storage,                           &
                                        Last_Column_Storage






USE Additional_Functions_Module, &
                                ONLY :  Lagrange_Poly,                                  &
                                        Spherical_Harmonic,                             &
                                        Map_To_X_Space, Map_From_X_Space,               &
                                        MVMULT_FULL

USE Functions_Quadrature, &
                                ONLY :  Initialize_LGL_Quadrature,                      &
                                        Initialize_LGL_Quadrature_Locations,            &
                                        Initialize_LG_Quadrature_Locations




USE Mesh_Module, &
                                ONLY :  Allocate_Mesh,                                  &
                                        Deallocate_Mesh,                                &
                                        Generate_Defined_Mesh


USE Source_Vector_Module, &
                                ONLY :  Allocate_Source_Vector,                         &
                                        Deallocate_Source_Vector,                       &
                                        Generate_Source_Vector


USE Stiffness_Matrix_Module, &
                                ONLY :  Allocate_Stiffness_Matrix,                      &
                                        Deallocate_Stiffness_Matrix,                    &
                                        Initialize_Stiffness_Matrix_Values,             &
                                        Generate_Stiffness_Matrix


USE Coefficient_Vector_Module, &
                                ONLY :  Allocate_Coefficient_Vector,                    &
                                        Deallocate_Coefficient_Vector,                  &
                                        Calculate_Coefficient_Vector


USE Timers_Module, &
                                ONLY :  Init_Timers,                    &
                                        Finalize_Timers,                &
                                        TimerStart,                     &
                                        TimerStop,                      &
                                        Timer_Matrix_Construction,      &
                                        Timer_SrcVec_Construction,      &
                                        Timer_LinSlv_Total


USE IO_Print_Setup_Module, &
                                ONLY :  Print_Setup

IMPLICIT NONE


CONTAINS


!+101+####################################################################################!
!                                                                                           !
!       Poseidon_Initialize                                                                 !
!                                                                                           !
!===========================================================================================!
!                                                                                           !
!   Sets code parameters, allocates space, and initializes functions and variables needed   !
!   to run Poseidon.                                                                        !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       Input Variables     :                                                               !
!                                                                                           !
!       FEM_Degree_Input        -       Integer, Order of Finite Element Method solver to   !
!                                               be performed.                               !
!                                                                                           !
!       L_Limit_Input           -       Integer, Limit on the l value of the Spherical      !
!                                               Harmonic spectral decomposition.            !
!                                                                                           !
!       Inner_Radius            -       Real, Inner radius of the computational domain.     !
!                                                                                           !
!       Outer_Radius            -       Real, Outer radius of the computational domain.     !
!                                                                                           !
!       R_Elements_Input        -       Integer, Number of radial elements.                 !
!                                                                                           !
!       T_Elements_Input        -       Integer, Number of theta elements.                  !
!                                                                                           !
!       P_Elements_Input        -       Integer, Number of phi elements.                    !
!                                                                                           !
!       Input_Delta_R_Vector    -       Optional Real Vector, Dimension(1:R_Elements_Input) !
!                                       Each value corresponds to the radial length of the  !
!                                       each radial shell of elements starting from         !
!                                       Inner_Radius.                                       !
!                                                                                           !
!       Input_Delta_T_Vector    -       Optional Real Vector, Dimension(1:T_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       theta wedge of elements starting from 0.            !
!                                                                                           !
!       Input_Delta_P_Vector    -       Optional Real Vector, Dimension(1:P_Elements_Input) !
!                                       Each value correspons to the angular width of each  !
!                                       phi wedge of elements starting from 0.              !
!                                                                                           !
 !#########################################################################################!

SUBROUTINE Poseidon_Initialize(     FEM_Degree_Input, L_Limit_Input, Inner_Radius, Outer_Radius,            &
                                    R_Elements_Input, T_Elements_Input, P_Elements_Input,                   &
                                    Num_R_Quad_Input, Num_T_Quad_Input, Num_P_Quad_Input,                   &
                                    Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector            )



                                         !                          !
                                        !!      Input Variables     !!
                                         !

INTEGER, INTENT(IN)                                                             ::  FEM_Degree_Input,       &
                                                                                    L_Limit_Input,          &
                                                                                    R_Elements_Input,       &
                                                                                    T_Elements_Input,       &
                                                                                    P_Elements_Input,       &
                                                                                    Num_R_Quad_Input,       &
                                                                                    Num_T_Quad_Input,       &
                                                                                    Num_P_Quad_Input

REAL(KIND = idp), INTENT(IN)                                                    ::  Inner_Radius,           &
                                                                                    Outer_Radius



REAL(KIND = idp), DIMENSION(1:R_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector





                                         !                              !
                                        !!     Subroutine Variables     !!
                                         !                              !

REAL(KIND = idp), DIMENSION(1:R_Elements_Input)                                 ::  Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:T_Elements_Input)                                 ::  Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:P_Elements_Input)                                 ::  Delta_P_Vector



CALL Init_Timers()



 !                                          !
!!  Set Global Variables to Input Values    !!
 !                                          !
DEGREE = FEM_Degree_Input
L_LIMIT = L_Limit_Input
LM_Length = (L_Limit+1)*(L_Limit+1)

R_INNER = Inner_Radius
R_OUTER = Outer_Radius


NUM_R_ELEMENTS = R_Elements_Input
NUM_T_ELEMENTS = T_Elements_Input
NUM_P_ELEMENTS = P_Elements_Input

NUM_R_NODES = DEGREE*NUM_R_ELEMENTS + 1




Source_Degrees(1) = Num_R_Quad_Input
Source_Degrees(2) = Num_T_Quad_Input
Source_Degrees(3) = Num_P_Quad_Input


Num_Quad_Dof = Source_Degrees(1)        &
             * Source_Degrees(2)        &
             * Source_Degrees(3)



                                 !                                      !
                                !!      Allocate Space for Poseidon     !!
                                 !                                      !




!!!!  Allocate Data Space !!!!
CALL Allocate_Mesh()
CALL Allocate_Source_Vector()
CALL Allocate_Coefficient_Vector()
Call Allocate_Stiffness_Matrix()


ALLOCATE(Source_Term_Coefficients(0:NUM_R_ELEMENTS-1,0:NUM_T_ELEMENTS-1,0:NUM_P_ELEMENTS-1,          &
                                  1:Source_Degrees(1),1:Source_Degrees(2),1:Source_Degrees(3) )     )

ALLOCATE( Source_Terms(1:Num_Quad_Dof, 0:NUM_R_ELEMENTS-1,0:NUM_T_ELEMENTS-1,0:NUM_P_ELEMENTS-1))







                                         !                                                      !
                                        !!      Initialize Stiffness Matrix Reusable Values     !!
                                         !                                                      !



!!                  Initalize STF_MAT Integrals                     !!
!!                                                                  !!
!!  * Only Needs Redoing if Num_R_Elements and/or Degree Change *   !!
!!                                                                  !!
CALL TimerStart( Timer_Matrix_Construction )
CALL Initialize_Stiffness_Matrix_Values()
CALL TimerStop( Timer_Matrix_Construction )












                                         !                                          !
                                        !!      Set Initial Mesh (Optional)         !!
                                         !                                          !


!!                                                                                              !!
!!      The Delta_*_Vector input are optional.  If they are provided Poseidon can go ahead      !!
!!      construct the data space meshes here.  This option is provided with the intent of       !!
!!      allowing static mesh users the chance to initialize the mesh here and therefore not     !!
!!      have to make a second call to "Poseidon_Set_Mesh".                                      !!
!!                                                                                              !!




  !                                                                                             !
 !!     If Input_Delta_R_Vector is present then user wishes to define an initial radial mesh.   !!
  !                                                                                             !
IF (    PRESENT(Input_Delta_R_Vector)   ) THEN



    CALL Generate_Defined_Mesh(NUM_R_ELEMENTS, R_INNER, Input_Delta_R_Vector, rlocs)
    RADIAL_MESH_SET_FLAG = .TRUE.



    !!              Generate Initial Stiffness Matrix                   !!
    !!                                                                  !!
    !!    Needs to be called everytime the radial mesh locations are    !!
    !!  set/altered, therefore is automatically called everytime the    !!
    !!  subroutine, Poseidon_Set_Mesh, is called to alter the radial    !!
    !!  mesh after this initial call.                                   !!
    !!                                                                  !!
    CALL Generate_Stiffness_Matrix()


ELSE IF ( ( RADIAL_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_R_Vector) ) ) THEN


    Delta_R_Vector = (R_OUTER - R_INNER) / REAL(Num_R_Elements )

    CALL Generate_Defined_Mesh(NUM_R_ELEMENTS, R_INNER, Delta_R_Vector, rlocs)
    RADIAL_MESH_SET_FLAG = .TRUE.

END IF





  !                                                                                         !
 !!   If Input_Delta_T_Vector is present then user wishes to define a non-uniform theta     !!
!!!   mesh. If this vector is not present then the code generates a Delta_T_Vector that     !!!
 !!   will produce a uniform mesh on the radial domain given by [0,Pi].                     !!
  !                                                                                         !
IF (    PRESENT(Input_Delta_T_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, Input_Delta_T_Vector, tlocs)
    THETA_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( THETA_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_T_Vector) ) ) THEN

    Delta_T_Vector = (pi) / REAL(Num_T_Elements )


    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, Delta_T_Vector, tlocs)
    THETA_MESH_SET_FLAG = .TRUE.



END IF






  !                                                                                         !
 !!   If Input_Delta_P_Vector is present then user wishes to define a non-uniform phi       !!
!!!   mesh. If this vector is not present then the code generates a Delta_P_Vector that     !!!
 !!   will produce a uniform mesh on the radial domain given by [0,2*Pi].                   !!
  !                                                                                         !
IF (    PRESENT(Input_Delta_P_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, Input_Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( PHI_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_P_Vector) ) ) THEN

    Delta_P_Vector = (2*pi) / REAL(Num_P_Elements )

    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.

END IF



CALL Print_Setup()







END SUBROUTINE Poseidon_Initialize


END MODULE Poseidon_Initialize_Module
