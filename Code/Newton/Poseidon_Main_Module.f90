   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Poseidon_Main_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +i02+   Poseidon_Newtonian_Potential_Output                                 !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize                                                 !##!
!##!    +102+   Poseidon_Run                                                        !##!
!##!    +103+   Poseidon_Close                                                      !##!
!##!    +104+   Poseidon_Set_Mesh                                                   !##!
!##!    +105+   Poseidon_Set_Boundary_Condtion                                      !##!
!##!                                                                                !##
!##!                                                                                !##!
!##!    +301+   Poseidon_1D_Newtonian_Potential_Output                              !##!
!##!    +302+   Poseidon_2D_Newtonian_Potential_Output                              !##!
!##!    +303+   Poseidon_3D_Newtonian_Potential_Output                              !##!
!##!    +304+   Poseidon_Newtonian_Center_And_Radial_Faces_Output                   !##!
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
            ONLY :  idp, pi



USE Global_Variables_And_Parameters, &
                                ONLY :  DEGREE, L_LIMIT,                                &
                                        R_INNER, R_OUTER,                               &
                                        NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS, &
                                        NUM_R_NODES,                                    &
                                        rlocs, tlocs, plocs,                            &
                                        Coefficient_Vector,                             &
                                        Source_Term_Coefficients,                       &
                                        Source_Terms,                                    &
                                        Source_Degrees,                                 &
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
                                        Timer_Core_SourceVector,        &
                                        Timer_Core_LinearSolve

USE IO_Print_Setup_Module, &
                                ONLY :  Print_Setup

IMPLICIT NONE

                !*I*============================================!
                !                                               !
                !                   Interfaces                  !
                !                                               !
                !===============================================!



 !+i02+################################################################################!
!                                                                                       !
!   Poseidon_Newtonian_Potential_Output -   Interface allowing for simplified solution  !
!                                       retrival depending on dimensionality of source. !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Acceptable Call Forms (See Associated Routine for description of variables)         !
!                                                                                       !
!   1 Dimension :   Poseidon_1D_Newtonian_Potential_Output                              !
!                                                                                       !
!           Left_Limit      -       Integer                                             !
!           Right_Limit     -       Integer                                             !
!           Num_Nodes       -       Integer                                             !
!           Output_R_Quad   -       Real Vector, Dimension(1:Num_Nodes)                 !
!           Potential       -       Real Vector, Dimension( 1:Num_Nodes,        &       !
!                                                           1:Num_R_Elements,   &       !
!                                                           1,                  &       !
!                                                           1                   )       !
!                                                                                       !
!   2 Dimensions :  Poseidon_3D_Newtonian_Potential_Output                              !
!                                                                                       !
!           Left_Limit      -       Integer                                             !
!           Right_Limit     -       Integer                                             !
!           Num_Nodes       -       Integer Vector, Dimension(1:2)                      !
!           Output_R_Quad   -       Real Vector, Dimension(1:Num_Nodes(1))              !
!           Output_T_Quad   -       Real Vector, Dimension(1:Num_Nodes(2))              !
!           Potential       -       Real Vector, Dimension( 1:Num_Nodes,        &       !
!                                                           1:Num_R_Elements,   &       !
!                                                           1:Num_T_Elements,   &       !
!                                                           1                   )       !
!                                                                                       !
!   3 Dimensions :  Poseidon_3D_Newtonian_Potential_Output                              !
!                                                                                       !
!           Left_Limit      -       Integer                                             !
!           Right_Limit     -       Integer                                             !
!           Num_Nodes       -       Integer Vector, Dimension(1:3)                      !
!           Output_R_Quad   -       Real Vector, Dimension(1:Num_Nodes(1))              !
!           Output_T_Quad   -       Real Vector, Dimension(1:Num_Nodes(2))              !
!           Output_P_Quad   -       Real Vector, Dimension(1:Num_Nodes(3))              !
!           Potential       -       Real Vector, Dimension( 1:Num_Nodes,        &       !
!                                                           1:Num_R_Elements,   &       !
!                                                           1:Num_T_Elements,   &       !
!                                                           1:Num_P_Elements    )       !
!                                                                                       !
 !#####################################################################################!
INTERFACE Poseidon_Newtonian_Potential_Output

    PROCEDURE   Poseidon_1D_Newtonian_Potential_Output, &
                Poseidon_2D_Newtonian_Potential_Output, &
                Poseidon_3D_Newtonian_Potential_Output

END INTERFACE


                    !*F&S*==========================================!
                    !                                               !
                    !           Functions & Subroutines             !
                    !                                               !
                    !===============================================!
CONTAINS
!+102+#######################################################################################!
!                                                                                           !
!       Poseidon_Run - Calculates Source Vector and Solves for solution coefficients        !
!                                                                                           !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Run(Time_Teller)

REAL(KIND = idp),OPTIONAL, INTENT(INOUT)            :: Time_Teller


LOGICAL                                             :: Readiness_Flag
REAL(KIND = idp)                                    :: Time_Keeper, Time_Keepers_Friend

CALL Poseidon_Readiness_Check(Readiness_Flag)


IF ( Readiness_Flag ) THEN


        !!! Generate Src Vector !!!
        CALL TimerStart( Timer_Core_SourceVector )
        CALL Generate_Source_Vector()
        CALL TimerStop( Timer_Core_SourceVector )



        !!! Calculate Solution Coefficients !!!
        CALL TimerStart( Timer_Core_LinearSolve )
        CALL Calculate_Coefficient_Vector()
        CALL TimerStop( Timer_Core_LinearSolve )



ELSE



    PRINT*, "ERROR IN POSEIDON : There was an error in setting up Poseidon, therefore it did not run."



END IF





END SUBROUTINE Poseidon_Run

























!+103+######################################################################################!
!                                                                                           !
!       Poseidon_Close - Deallocate all Poseidon variables.                                 !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Close()



!!!!  Deallocate Data Space !!!!
CALL Deallocate_Mesh
CALL Deallocate_Source_Vector()
CALL Deallocate_Coefficient_Vector()
CALL Deallocate_Stiffness_Matrix()


DEALLOCATE(Source_Term_Coefficients)

IF (ALLOCATED(First_Column_Storage) ) THEN

    DEALLOCATE(First_Column_Storage)

END IF


IF (ALLOCATED(Last_Column_Storage) ) THEN

    DEALLOCATE(Last_Column_Storage)

END IF



IF (    Test_Space_Allocated_Flag .EQV. .TRUE.  ) THEN

    DEALLOCATE(Test_Source_Input)

END IF



Stiffness_Matrix_Initialized_Flag = .FALSE.
Test_Space_Allocated_Flag = .FALSE.
Test_Run_Flag = .FALSE.


RADIAL_MESH_SET_FLAG = .FALSE.
THETA_MESH_SET_FLAG = .FALSE.
PHI_MESH_SET_FLAG = .FALSE.

INNER_BC_SET_FLAG = .FALSE.
OUTER_BC_SET_FLAG = .FALSE.


CALL Finalize_Timers()

END SUBROUTINE Poseidon_Close
















 !+104+####################################################################################!
!                                                                                           !
!      Poseidon_Set_Mesh                                                                    !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   This subroutine sets the locations of the element end locations using vectors giving    !
!   the length of elements in each dimension.                                               !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :      * = R, T, P                                                  !
!                                                                                           !
!   Input_Delta_*_Vector    -   Optional, Real Vector, Dimension(1:NUM_*_ELEMENTS)          !
!                               Gives the length of elements in the * dimension.            !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Set_Mesh(Input_Delta_R_Vector, Input_Delta_T_Vector, Input_Delta_P_Vector)



REAL(KIND = idp), DIMENSION(1:NUM_R_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:NUM_T_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:NUM_P_ELEMENTS),    OPTIONAL,   INTENT(IN)      ::  Input_Delta_P_Vector

REAL(KIND = idp), DIMENSION(1:NUM_R_ELEMENTS)                                  ::  TMP_Delta_R_Vector
REAL(KIND = idp), DIMENSION(1:NUM_T_ELEMENTS)                                  ::  TMP_Delta_T_Vector
REAL(KIND = idp), DIMENSION(1:NUM_P_ELEMENTS)                                  ::  TMP_Delta_P_Vector



  !                                                                                                 !
 !!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
  !                                                                                                 !
IF (    PRESENT(Input_Delta_R_Vector)   ) THEN




    CALL Generate_Defined_Mesh(NUM_R_ELEMENTS, R_INNER, Input_Delta_R_Vector, rlocs)
    RADIAL_MESH_SET_FLAG = .TRUE.

     !                                                                              !
    !!  If the stiffness matrix has been initialized then this call is resetting    !!
    !!  the radial mesh, which means the stiffness matrix needs to be regenerated.  !!
     !                                                                              !
    IF ( Stiffness_Matrix_Initialized_Flag .EQV. .TRUE. ) THEN

        CALL Generate_Stiffness_Matrix()

    END IF




END IF





 !                                                                                                 !
!!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
 !                                                                                                 !
IF (    PRESENT(Input_Delta_T_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, Input_Delta_T_Vector, tlocs)
    THETA_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( THETA_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_T_Vector) ) ) THEN


    TMP_Delta_R_Vector = (pi)/ REAL(NUM_T_ELEMENTS)

    CALL Generate_Defined_Mesh(NUM_T_ELEMENTS, 0.0_idp, (/ pi /), tlocs)
    THETA_MESH_SET_FLAG = .TRUE.



END IF





 !                                                                                                 !
!!     If Input_Delta_R_Vector is present then user wishes to define a non-uniform radial mesh.    !!
 !                                                                                                 !
IF (    PRESENT(Input_Delta_P_Vector)   ) THEN


    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, Input_Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.


ELSE IF ( ( PHI_MESH_SET_FLAG .EQV. .FALSE. ) .AND. ( .NOT. PRESENT(Input_Delta_P_Vector) ) ) THEN

    TMP_Delta_P_Vector = (2*pi) / REAL(Num_P_Elements )


    CALL Generate_Defined_Mesh(NUM_P_ELEMENTS, 0.0_idp, TMP_Delta_P_Vector, plocs)
    PHI_MESH_SET_FLAG = .TRUE.

END IF









END SUBROUTINE Poseidon_Set_Mesh













 !+105+####################################################################################!
!                                                                                           !
!      Poseidon_Set_Boundary_Condtion                                                       !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Set boundary condition values for the system.  This needs to be done once, but can be   !
!   repeatedly used if the boundary conditions are changed.                                 !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables     :                                                                   !
!                                                                                           !
!           BC_Location_Input       -   "I" for Inner Boundary                              !
!                                       "O" for Outer Boundary                              !
!                                                                                           !
!           BC_Type_Input           -   "D" for Dirichlet Boundary Condition                !
!                                       "N" for Neumann Boundary Condition                  !
!                                                                                           !
!           BC_Value_Input          -   For a Dirichlet Boundary Condition specify the      !
!                                           Newtonian potential at the boundary.            !
!                                       For a Neumann Boundary Condition specify the        !
!                                           radial derivative of the potential at the       !
!                                           boundary.                                       !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Set_Uniform_Boundary_Condition(BC_Location_Input, BC_Type_Input, BC_Value_Input)



CHARACTER(LEN = *),         INTENT(IN)                  ::  BC_Location_Input
CHARACTER(LEN = *),         INTENT(IN)                  ::  BC_Type_Input

REAL(KIND = idp),           INTENT(IN)                  ::  BC_Value_Input



!                               !
!   Inner Boundary Condition    !
!                               !
IF (    BC_Location_Input == "I"    ) THEN




    !                                   !
    !   Dirichlet Boundary Condition    !
    !                                   !
    IF ( (  BC_Type_Input == "D"  )  ) THEN


        INNER_DIR_BC_INPUT = BC_Value_Input
        INNER_BC_TYPE = "DRCH"




    !                                   !
    !    Neumann Boundary Condition     !
    !                                   !
    ELSE IF ( (  BC_Type_Input == "N"  )  ) THEN


        INNER_NEU_BC_INPUT = BC_Value_Input
        INNER_BC_TYPE = "NEUM"


    END IF


    INNER_BC_SET_FLAG = .TRUE.
    INNER_UNIFORM_DIR_BC_FLAG = .TRUE.






!                               !
!   Outer Boundary Condition    !
!                               !
ELSE IF (   BC_Location_Input == "O"    )  THEN




    !                                   !
    !   Dirichlet Boundary Condition    !
    !                                   !
    IF ( (  BC_Type_Input == "D"  ) .OR. (  BC_Type_Input == "D"  ) ) THEN


        OUTER_DIR_BC_INPUT = BC_Value_Input
        OUTER_BC_TYPE = "DRCH"


    !                                   !
    !   Dirichlet Boundary Condition    !
    !                                   !
    ELSE IF ( (  BC_Type_Input == "N"  ) .OR. (  BC_Type_Input == "N"  ) ) THEN


        OUTER_NEU_BC_INPUT = BC_Value_Input
        OUTER_BC_TYPE = "NEUM"


    END IF


    OUTER_BC_SET_FLAG = .TRUE.
    OUTER_UNIFORM_DIR_BC_FLAG = .TRUE.


END IF








END SUBROUTINE Poseidon_Set_Uniform_Boundary_Condition







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







 !+301+####################################################################################!
!                                                                                           !
!                           Poseidon_1D_Newtonian_Potential_Output                          !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides an simplified output method for the result of Poseidon's   !
!   calulations. It takes in a list of locations and returns the potential value at that    !
!   location by using the expansion coefficients to reconstruct the solution.               !
!                                                                                           !
!       The 1-Dimensional version of this subroutine takes in variables describing the      !
!   locations where the solution is desired, and sends them along with filler variables to  !
!   the subroutine, Poseidon_3D_Newtonian_Potential_Output (+303+). The filler variables    !
!   instruct that subroutine to act only in one dimension.                                  !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables                                                                         !
!                                                                                           !
!   Left_Limit,         -   Integer, Left limit of the output quadrature space.             !
!                                                                                           !
!   Right_Limit         -   Integer, Right limit of the output quadrature space.            !
!                                                                                           !
!   Num_Nodes           -   Integer, Number of radial nodes                                 !
!                                                                                           !
!   Output_R_Quad       -   Real Vector, Locations of radial nodes in the output quadrature !
!                               space.                                                      !
!                                                                                           !
!                                                                                           !
!                                                                                           !
!   Output Variables                                                                        !
!                                                                                           !
!   Potential           -   Real Vector, Potential value at each location in each element   !
!                                                                                           !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_1D_Newtonian_Potential_Output(              Left_Limit,             &
                                                                Right_Limit,            &
                                                                Num_Nodes,              &
                                                                Output_R_Quad,          &
                                                                Potential                   )





                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp),                               INTENT(IN)      ::  Left_Limit,             &
                                                                    Right_Limit


INTEGER,                                        INTENT(IN)      ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes),       INTENT(IN)      ::  Output_R_Quad


REAL(KIND = idp), DIMENSION(1:Num_Nodes,                        &
                            1:NUM_R_ELEMENTS,                   &
                            1:NUM_T_ELEMENTS,                   &
                            1:NUM_P_ELEMENTS),  INTENT(INOUT)   ::  Potential





                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !

INTEGER, DIMENSION(1:3)                                         ::  Num_Nodes_Vectorized

REAL(KIND = idp), DIMENSION(1:1)                                ::  Output_T_Quad
REAL(KIND = idp), DIMENSION(1:1)                                ::  Output_P_Quad





 !                                                                   !
!!   The 3D subroutine requires a length 3 vector for Num_Nodes.     !!
 !                                                                   !
Num_Nodes_Vectorized(1) = Num_Nodes
Num_Nodes_Vectorized(2) = 1
Num_Nodes_Vectorized(3) = 1




   !                                                                     !
  !!   Sets the location of the single node in the Theta and Phi         !!
 !!!   dimensions to 1.  This choice is arbitrary as this dimension is   !!!
!!!!   being dismissed, which is achieved by using the lowert order      !!!!
 !!!   Lagrange Polynomial, L(x) = 1, which is independent of input      !!!
  !!   location.                                                         !!
   !                                                                     !
Output_T_Quad = 1.0_idp
Output_P_Quad = 1.0_idp









CALL Poseidon_3D_Newtonian_Potential_Output(Left_Limit, Right_Limit, Num_Nodes_Vectorized,          &
                                            Output_R_Quad, Output_T_Quad, Output_P_Quad, Potential      )








END SUBROUTINE Poseidon_1D_Newtonian_Potential_Output









 !+302+####################################################################################!
!                                                                                           !
!                           Poseidon_2D_Newtonian_Potential_Output                          !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides an simplified output method for the result of Poseidon's   !
!   calulations. It takes in a list of locations and returns the potential value at that    !
!   location by using the expansion coefficients to reconstruct the solution.               !
!                                                                                           !
!       The 2-Dimensional version of this subroutine takes in variables describing the      !
!   locations where the solution is desired, and sends them along with filler variables to  !
!   the subroutine, Poseidon_3D_Newtonian_Potential_Output (+303+). The filler variables    !
!   instruct that subroutine to act only in two dimensions.                                 !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables                                                                         !
!                                                                                           !
!   Left_Limit,         -   Integer, Left limit of the output quadrature space.             !
!                                                                                           !
!   Right_Limit         -   Integer, Right limit of the output quadrature space.            !
!                                                                                           !
!   Num_Nodes           -   Integer Vector, Number of nodes per dimension                   !
!                                                                                           !
!   Output_*_Quad       -   Real Vectors, Locations of nodes in each dimension of the       !
!                               output quadrature space.                                    !
!                                                                                           !
!                                                                                           !
!                                                                                           !
!   Output Variables                                                                        !
!                                                                                           !
!   Potential           -   Real Vector, Potential value at each location in each element   !
!                                                                                           !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_2D_Newtonian_Potential_Output(                  Left_Limit,             &
                                                                    Right_Limit,            &
                                                                    Num_Nodes,              &
                                                                    Output_R_Quad,          &
                                                                    Output_T_Quad,          &
                                                                    Potential                   )


                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp),                               INTENT(IN)          ::  Left_Limit,         &
                                                                        Right_Limit


INTEGER, DIMENSION(1:2),                        INTENT(IN)          ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)),    INTENT(IN)          ::  Output_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2)),    INTENT(IN)          ::  Output_T_Quad  !??? Right Dimension ???!


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2),            &
                            1:NUM_R_ELEMENTS,                       &
                            1:NUM_T_ELEMENTS,                       &
                            1:NUM_P_ELEMENTS),  INTENT(INOUT)       ::  Potential







                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !

INTEGER, DIMENSION(1:3)                                         ::  Num_Nodes_Vectorized


REAL(KIND = idp), DIMENSION(1:1)                                ::  Output_P_Quad






 !                                                                   !
!!   The 3D subroutine requires a length 3 vector for Num_Nodes.     !!
 !                                                                   !
Num_Nodes_Vectorized(1) = Num_Nodes(1)
Num_Nodes_Vectorized(2) = Num_Nodes(2)
Num_Nodes_Vectorized(3) = 1



  !                                                                     !
 !!   Sets the location of the single node in the Phi dimension to      !!
!!!   1.  This choice is arbitrary as this dimension is being           !!!
!!!   dismissed, which is achieved by using the lowert order Lagrange   !!!
 !!   Polynomial, L(x) = 1, which is independent of input location.     !!
  !                                                                     !
Output_P_Quad = 1.0_idp







CALL Poseidon_3D_Newtonian_Potential_Output(Left_Limit, Right_Limit, Num_Nodes_Vectorized,          &
                                            Output_R_Quad, Output_T_Quad, Output_P_Quad, Potential      )











END SUBROUTINE Poseidon_2D_Newtonian_Potential_Output














 !+303+####################################################################################!
!                                                                                           !
!                           Poseidon_3D_Newtonian_Potential_Output                          !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides an output method for the result of Poseidon's calulations. !
!   It takes in a list of locations and returns the potential value at that location by     !
!   using the expansion coefficients to reconstruct the solution.                           !
!                                                                                           !
!       This subroutine is capable of handling 1- and 2-Dimensional outputs if the input    !
!   is properly established. To simplify the input of lower dimensional systems, the        !
!   subroutines, Poseidon_1D_Newtonain_Potential_Output (+301+), and                        !
!   Poseidon_2D_Newtonain_Potential_Output (+302+) have been created.  These subroutines    !
!   require less input, and generate filler variables that limit this subroutine actions    !
!   to the appropriote dimension.                                                           !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables                                                                         !
!                                                                                           !
!   Left_Limit,         -   Integer, Left limit of the output quadrature space.             !
!                                                                                           !
!   Right_Limit         -   Integer, Right limit of the output quadrature space.            !
!                                                                                           !
!   Num_Nodes           -   Integer Vector, Number of nodes per dimension                   !
!                                                                                           !
!   Output_*_Quad       -   Real Vectors, Locations of nodes in each dimension of the       !
!                               output quadrature space.                                    !
!                                                                                           !
!                                                                                           !
!                                                                                           !
!   Output Variables                                                                        !
!                                                                                           !
!   Potential           -   Real Vector, Potential value at each location in each element   !
!                                                                                           !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_3D_Newtonian_Potential_Output(              Left_Limit,                 &
                                                                Right_Limit,                &
                                                                Num_Nodes,                  &
                                                                Output_R_Quad,              &
                                                                Output_T_Quad,              &
                                                                Output_P_Quad,              &
                                                                Potential                       )


                             !                  !
                            !!  Input Variables !!
                             !                  !
REAL(KIND = idp), INTENT(IN)                                            ::  Left_Limit,    &
                                                                            Right_Limit


INTEGER, DIMENSION(1:3), INTENT(IN)                                     ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)),    INTENT(IN)              ::  Output_R_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2)),    INTENT(IN)              ::  Output_T_Quad
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3)),    INTENT(IN)              ::  Output_P_Quad



REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:NUM_R_ELEMENTS,                           &
                            1:NUM_T_ELEMENTS,                           &
                            1:NUM_P_ELEMENTS),  INTENT(OUT)             ::  Potential










                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !
INTEGER                                                     ::  re, te, pe,                     &
                                                                Output_R, Output_T, Output_P,   &
                                                                l, m, d,                        &
                                                                Output_Here,                    &
                                                                Num_Output_DOF


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_X_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2))                 ::  Output_T_X_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3))                 ::  Output_P_X_Locations

REAL(KIND = idp), DIMENSION(1:Num_Nodes(2))                 ::  Output_T_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3))                 ::  Output_P_Locations


REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Local_Locations
REAL(KIND = idp), DIMENSION(0:DEGREE)                       ::  Lagrange_Poly_Values



COMPLEX(KIND = idp)                                         ::  Tmp_Sphere_Harmonic_Value,      &
                                                                Tmp_Potential_Value

                             !                    !
                            !!  Global Variables  !!
                             !                    !

!INTEGER,                                                   ::  NUM_R_ELEMENTS,             &
!                                                               NUM_T_ELEMENTS,             &
!                                                               NUM_P_ELEMENTS









                          !                                                 !
                         !!                                                 !!
                        !!!               Calculate Local DOF               !!!
                         !!                                                 !!
                          !                                                 !

Num_Output_DOF = Num_Nodes(1) * Num_Nodes(2) * Num_Nodes(3)







                          !                                                 !
                         !!                                                 !!
                        !!!          Initialize Local Quadratures           !!!
                         !!                                                 !!
                          !                                                 !

Local_Locations = Initialize_LGL_Quadrature_Locations(DEGREE)





                          !                                                 !
                         !!                                                 !!
                        !!!     Map Output Locations to [ -1, 1 ] Space      !!!
                         !!                                                 !!
                          !                                                 !

Output_R_X_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Output_R_Quad)
Output_T_X_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Output_T_Quad)
Output_P_X_Locations = Map_To_X_Space(Left_Limit, Right_Limit, Output_P_Quad)










DO re = 0,NUM_R_ELEMENTS-1


     !                                                                                      !
    !!   No need to map radial locations, as Lagrange Polynomial works in [-1, 1] space     !!
     !                                                                                      !


    DO te = 0,NUM_T_ELEMENTS-1



         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        Output_T_Locations = Map_From_X_Space(tlocs(te), tlocs(te + 1), Output_T_X_Locations)


        DO pe = 0,NUM_P_ELEMENTS-1



             !                                                          !
            !!    Map Phi Locations from [-1,1] space to real space.    !!
             !                                                          !
            Output_P_Locations = Map_From_X_Space(plocs(pe), plocs(pe + 1), Output_P_X_Locations)



             !                                          !
            !!   Set/Reset Output Vector Location       !!
             !                                          !
            Output_Here = 1


            DO Output_P = 1,Num_Nodes(3)
            DO Output_T = 1,Num_Nodes(2)
            DO Output_R = 1,Num_Nodes(1)


                 !                                                                      !
                !!   Calculate Lagrange Polynomial Values for current output location   !!
                 !                                                                      !
                Lagrange_Poly_Values = Lagrange_Poly(Output_R_X_Locations(Output_R), DEGREE, Local_Locations)


                 !                                                  !
                !!   Set/Reset temporary value holder to zero.      !!
                 !                                                  !
                Tmp_Potential_Value = 0.0_idp





                DO l = 0,L_Limit
                DO m = -l,l

                     !                                                                          !
                    !!  Precompute and store the value of the spherical harmonic function for   !!
                    !!  the current location and choice of l, and m.                            !!
                     !                                                                          !
                    Tmp_Sphere_Harmonic_Value = Spherical_Harmonic( l,                                  &
                                                                m,                                  &
                                                                Output_T_Locations(Output_T),       &
                                                                Output_P_Locations(Output_P)        )



                    DO d = 0,DEGREE



                         !                                                                  !
                        !!   Reconstruct the potential using the expansion coefficients     !!
                         !                                                                  !
                        Tmp_Potential_Value = Tmp_Potential_Value + Coefficient_Vector(re*DEGREE + d, m, l)     &
                                                                  * Tmp_Sphere_Harmonic_Value                   &
                                                                  * Lagrange_Poly_Values(d)



                    END DO  !   d Loop
                END DO  !   m Loop
                END DO  !   l Loop


                 !                                                      !
                !!   Output only real part of solution to out vector    !!
                 !                                                      !
                Potential(Output_Here, re+1, te+1, pe+1) = REAL(Tmp_Potential_Value,idp)


                 !                                      !
                !!   Increment output vector location   !!
                 !                                      !
                Output_Here = Output_Here + 1





            END DO  !   Output_R Loop
            END DO  !   Output_T Loop
            END DO  !   Output_P Loop







        END DO  ! pe Loop

    END DO  ! te Loop

END DO  !   re Loop






END SUBROUTINE Poseidon_3D_Newtonian_Potential_Output








 !+304+####################################################################################!
!                                                                                           !
!                           Poseidon_Newtonian_Center_And_Radial_Faces_Output               !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!       This subroutine provides an simplified output method for the result of Poseidon's   !
!   calulations. It takes in a list of locations and returns the potential value at that    !
!   location by using the expansion coefficients to reconstruct the solution.               !
!                                                                                           !
!       The 1-Dimensional version of this subroutine takes in variables describing the      !
!   locations where the solution is desired, and sends them along with filler variables to  !
!   the subroutine, Poseidon_3D_Newtonian_Potential_Output (+303+). The filler variables    !
!   instruct that subroutine to act only in one dimension.                                  !
!                                                                                           !
!-------------------------------------------------------------------------------------------!
!                                                                                           !
!   Input Variables                                                                         !
!                                                                                           !
!   Left_Limit,         -   Integer, Left limit of the output quadrature space.             !
!                                                                                           !
!   Right_Limit         -   Integer, Right limit of the output quadrature space.            !
!                                                                                           !
!   Num_Nodes           -   Integer, Number of radial nodes                                 !
!                                                                                           !
!   Output_R_Quad       -   Real Vector, Locations of radial nodes in the output quadrature !
!                               space.                                                      !
!                                                                                           !
!                                                                                           !
!                                                                                           !
!   Output Variables                                                                        !
!                                                                                           !
!   Potential           -   Real Vector, Potential value at each location in each element   !
!                                                                                           !
!                                                                                           !
 !#########################################################################################!
SUBROUTINE Poseidon_Newtonian_Center_And_Radial_Faces_Output(       Potential       )





                             !                  !
                            !!  Input Variables !!
                             !                  !



REAL(KIND = idp), DIMENSION(1:3,                                &
                            1:NUM_R_ELEMENTS,                   &
                            1:NUM_T_ELEMENTS,                   &
                            1:NUM_P_ELEMENTS),  INTENT(INOUT)   ::  Potential





                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !

REAL(KIND = idp)                                                ::  Left_Limit,                 &
                                                                    Right_Limit


INTEGER, DIMENSION(1:3)                                         ::  Num_Nodes_Vectorized

REAL(KIND = idp), DIMENSION(1:3)                                ::  Output_R_Quad
REAL(KIND = idp), DIMENSION(1:1)                                ::  Output_T_Quad
REAL(KIND = idp), DIMENSION(1:1)                                ::  Output_P_Quad







 !                                                                   !
!!   The 3D subroutine requires a length 3 vector for Num_Nodes.     !!
 !                                                                   !
Num_Nodes_Vectorized(1) = 3     ! One node per face, and one for the center
Num_Nodes_Vectorized(2) = 1
Num_Nodes_Vectorized(3) = 1




  !                                                                         !
 !!     Sets the location of the single node in the Theta and Phi           !!
!!!     dimensions to 0.0 to select the angular center of the element.      !!!
!!!     The radial locations are chosen to select the two ends of the       !!!
 !!     the element (-1,1) and the center (0).                              !!
  !                                                                         !
Left_Limit = -1.0_idp
Right_Limit = 1.0_idp

Output_R_Quad = (/ Left_Limit, 0.0_idp, Right_Limit /)
Output_T_Quad = 0.0_idp
Output_P_Quad = 0.0_idp









CALL Poseidon_3D_Newtonian_Potential_Output(Left_Limit, Right_Limit, Num_Nodes_Vectorized,          &
                                            Output_R_Quad, Output_T_Quad, Output_P_Quad, Potential      )








END SUBROUTINE Poseidon_Newtonian_Center_And_Radial_Faces_Output

















!+401+######################################################################################!
!                                                                                           !
!       Poseidon_Readiness_Check                                                            !
!                                                                                           !
!###########################################################################################!
SUBROUTINE Poseidon_Readiness_Check(Readiness_Flag)


LOGICAL,                                INTENT(INOUT)       :: Readiness_Flag
LOGICAL                                                     :: Error_Flag


Error_Flag = .FALSE.

                                    !               !
                                    !   Mesh Check  !
                                    !               !
IF ( .NOT. RADIAL_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The radial mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF

IF ( .NOT. THETA_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The theta mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF

IF ( .NOT. PHI_MESH_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The phi mesh was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF





                            !                               !
                            !   Boundary Conditions Check   !
                            !                               !
IF ( .NOT. INNER_BC_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : The inner boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF


IF ( .NOT. OUTER_BC_SET_FLAG ) THEN

    PRINT*,"ERROR IN POSEIDON : Ther outer boundary condition was not set before 'Poseidon_Run' was called."
    Error_Flag = .TRUE.

END IF













Readiness_Flag = .NOT. Error_Flag





END SUBROUTINE Poseidon_Readiness_Check





END MODULE Poseidon_Main_Module
