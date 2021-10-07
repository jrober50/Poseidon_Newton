   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
PROGRAM Poseidon_Newtonian                                                          !##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!    This program is an example of how to run the Newtonian side of Poseidon.    !##!
!##!                                                                                !##!
!##!
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
            ONLY :  idp, pi

USE COMMON_IO, &
            ONLY :  FILE_DATA_TYPE,         &
                    OPEN_FILE,              &
                    CLOSE_FILE



USE Global_Variables_And_Parameters,&
            ONLY :  R_INNER,                    &
                    R_OUTER,                    &
                    Analytic_Solution,          &
                    rlocs, plocs, tlocs,        &
                    Source_Term_Coefficients,   &
                    STF_ELEM_VAL, STF_COL_PTR, STF_ROW_IND,           &
                    Coefficient_Vector,         &
                    Source_Vector,              &
                    NUM_R_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    NUM_R_NODES,                &
                    Test_Source_Input




USE Poseidon_Main_Module, &
            ONLY :  Poseidon_Run,                           &
                    Poseidon_Close,                         &
                    Poseidon_Set_Uniform_Boundary_Condition

USE Poseidon_Source_Input_Module, &
            ONLY : Poseidon_Newtonian_Source_Input

USE Poseidon_Initialize_Module, &
            ONLY :  Poseidon_Initialize


USE Additional_Functions_Module, &
            ONLY :  Map_From_X_Space
USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature_Locations

USE Test_Functions_Module, &
            ONLY :  Poseidon_Initialize_Test_Problem


USE Timers_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop,                      &
                    Timer_Initialize_Test_Problem,  &
                    Timer_Core_SourceInput,         &
                    Timer_Core_PrintResults

USE IO_Print_Results_Module, &
            ONLY :  Print_Results


IMPLICIT NONE






                !*I*============================================!
                !                                               !
                !            Variable Initialization            !
                !                                               !
                !===============================================!

INTEGER                                                    ::   i, re, te, pe          ! DO Loop Counter Variables

INTEGER                                                     ::  FEM_Degree_Input,       &
                                                                L_Limit_Input

REAL(KIND = idp)                                            ::  Inner_Radius,           &
                                                                Outer_Radius


INTEGER                                                     ::  R_Elements_Input,      &
                                                                T_Elements_Input,      &
                                                                P_Elements_Input

INTEGER                                                     ::  Num_R_Quad_Input,      &
                                                                Num_T_Quad_Input,      &
                                                                Num_P_Quad_Input


REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_Delta_P_Vector




!                           !
!   Source Input Variables  !
!                           !


REAL(KIND = idp)                                            ::  Left_Limit,             &
                                                                Right_Limit


INTEGER, DIMENSION(1:3)                                     ::  Num_Input_Nodes

REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_R_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_T_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Input_P_Quad
REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Rho




!                               !
!   Boundary Value Variables    !
!                               !

REAL(KIND = idp)                                            :: Enclosed_Mass




!                       !
!   Output Variables    !
!                       !


REAL(KIND = idp)                                            ::  potential


REAL(KIND = idp)                                            ::  Output_Left_Limit,  &
                                                                Output_Right_Limit

INTEGER, DIMENSION(1:3)                                     ::  Num_Output_Nodes

REAL(KIND = idp), DIMENSION(:,:,:,:), ALLOCATABLE           ::  Potential_Output



REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_R_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_T_Quad
REAL(KIND = idp), DIMENSION(:), ALLOCATABLE                 ::  Output_P_Quad




!                                                       !
!   Other Variables Used to Initialize or Store Values  !
!                                                       !
INTEGER                                                     ::  Num_DOF, Num_Out_DOF





FEM_Degree_Input = 5
L_Limit_Input = 0


Inner_Radius = 0.0_idp
Outer_Radius = 1.0_idp


R_Elements_Input = 100
T_Elements_Input = 10                   !   For a 1-Dimensional Simulation Set equal to 1
P_Elements_Input = 10                   !   For a 1 or 2-Dimensional Simulations Set equal to 1


Num_R_Quad_Input = 5
Num_T_Quad_Input = 5                    !   For a 1-Dimensional Simulation Set equal to 1
Num_P_Quad_Input = 5                    !   For a 1 or 2-Dimensional Simulations Set equal to 1


ALLOCATE(Input_Delta_R_Vector(1:R_Elements_Input))
ALLOCATE(Input_Delta_T_Vector(1:T_Elements_Input))
ALLOCATE(Input_Delta_P_Vector(1:P_Elements_Input))


!                                               !
!   For this example, we build uniform meshes.  !
!                                               !
Input_Delta_R_Vector = (Outer_Radius - Inner_Radius)/REAL(R_Elements_Input)
Input_Delta_T_Vector = (pi)/REAL(T_Elements_Input)
Input_Delta_P_Vector = (2*pi)/REAL(P_Elements_Input)










Left_Limit = -0.5_idp
Right_Limit = 0.5_idp

!Num_Input_Nodes(1) = 5
!Num_Input_Nodes(2) = 1
!Num_Input_Nodes(3) = 1
Num_Input_Nodes(1) = Num_R_Quad_Input
Num_Input_Nodes(2) = Num_T_Quad_Input
Num_Input_Nodes(3) = Num_P_Quad_Input

Num_DOF = Num_Input_Nodes(1)*Num_Input_Nodes(2)*Num_Input_Nodes(3)



ALLOCATE(Input_R_Quad(1:Num_Input_Nodes(1)), Input_T_Quad(1:Num_Input_Nodes(2)), Input_P_Quad(1:Num_Input_Nodes(3)))
ALLOCATE(Rho(1:Num_DOF, 1:R_Elements_Input, 1:T_Elements_Input, 1:P_Elements_Input))



!                                                                       !
!   For this example the quadrature points in each dimension are the    !
!   Lobatto quadrature points.                                          !
!                                                                       !
Input_R_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(1))
Input_T_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(2))
Input_P_Quad = Initialize_LG_Quadrature_Locations(Num_Input_Nodes(3))




!                                                                       !
!   As the Lobatto points are defined in [-1,1] space, but we want to   !
!   provide source input on a [-.5, .5] space, we need to map the       !
!   the locations from one space to the other.                          !
!                                                                       !

!   Function - Map_From_X_Space - This Function maps Input_R_Quad (Real Double Value/Vector in [-1,1] space) to
!                                   a space defined by the limits (Real Double Values)
Input_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_R_Quad)
Input_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_T_Quad)
Input_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Input_P_Quad)









CALL Poseidon_Initialize(   FEM_Degree_Input,       &
                            L_Limit_Input,          &
                            Inner_Radius,           &
                            Outer_Radius,           &
                            R_Elements_Input,       &
                            T_Elements_Input,       &
                            P_Elements_Input,       &
                            Num_R_Quad_Input,       &
                            Num_T_Quad_Input,       &
                            Num_P_Quad_Input,       &
                            Input_Delta_R_Vector)!, &       ! The following two input vectors are optional.
!                           Input_Delta_T_Vector,   &       ! Does not need to be included for 1D problems
!                           Input_Delta_P_Vector         )  ! Does not need to be includes for 1D or 2D problems







!                                         !
!!                                       !!
!!!         Create Source Values        !!!
!!                                       !!
!                                         !


    !!!  Call To Run Test Problems !!!
CALL Poseidon_Initialize_Test_Problem(1, Num_Input_Nodes, Rho)












!                                         !
!!                                       !!
!!!       Set Boundary Conditions       !!!
!!                                       !!
!                                         !

Enclosed_Mass = (4.0_idp/3.0_idp)*pi!*0.5*0.5*0.5

CALL Poseidon_Set_Uniform_Boundary_Condition("I","N",0.0_idp)
CALL Poseidon_Set_Uniform_Boundary_Condition("O","D",-Enclosed_Mass/R_OUTER)










!                                         !
!!                                       !!
!!!        Input Source Values          !!!
!!                                       !!
!                                         !

CALL TimerStart( Timer_Core_SourceInput )

CALL Poseidon_Newtonian_Source_Input(       Left_Limit,             &
                                            Right_Limit,            &
                                            Num_Input_Nodes(1),           &
                                            Input_R_Quad,           &
                                            Rho                         )


CALL TimerStop( Timer_Core_SourceInput )






!                                         !
!!                                       !!
!!!             Run Poseidon            !!!
!!                                       !!
!                                         !
                                                                                          !
CALL Poseidon_Run()











!
! Output Results
!
Call TimerStart( Timer_Core_PrintResults )
IF ( 1 == 1 ) THEN

    Call Print_Results()

END IF
Call TimerStop( Timer_Core_PrintResults )






!
! Close Poseidon
!

CALL Poseidon_Close()









END PROGRAM Poseidon_Newtonian
