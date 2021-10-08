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
            ONLY :  Poseidon_Initialize,                    &
                    Poseidon_Run,                           &
                    Poseidon_Close,                         &
                    Poseidon_Newtonian_Source_Input,        &
                    Calculate_Potential_At_Location,        &
                    Poseidon_Newtonian_Potential_Output,    &
                    Poseidon_Set_Mesh,                      &
                    Poseidon_Newtonian_Center_And_Radial_Faces_Output,  &
                    Poseidon_Set_Uniform_Boundary_Condition




USE Additional_Functions_Module, &
            ONLY :  Initialize_LG_Quadrature_Locations,     &
                    Map_From_X_Space


USE Test_Functions_Module, &
            ONLY :  Poseidon_Initialize_Test_Problem



                    




IMPLICIT NONE






                !*I*============================================!
                !                                               !
                !            Variable Initialization            !
                !                                               !
                !===============================================!

INTEGER                                                    ::   i, re, te, pe          ! DO Loop Counter Variables




!                                       !
!   Poseidon Initialization Variables   !
!                                       !




INTEGER                                                     ::  FEM_Degree_Input,       &
                                                                L_Limit_Input

REAL(KIND = idp)                                            ::  Inner_Radius,           &
                                                                Outer_Radius


INTEGER                                                     ::  R_Elements_Input,      &
                                                                T_Elements_Input,      &
                                                                P_Elements_Input

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

INTEGER                                                     ::  NUM_SAMPLES
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

REAL(KIND = idp)                                            ::  deltar, r, theta, phi
INTEGER                                                     ::  Num_DOF, Num_Out_DOF
REAL(KIND = idp)                                            ::  TMP

















!!                                       !!
!!   Poseidon Initialization Variables   !!
!!                                       !!




!   Set the Degree of the Finite Element Expansion !
FEM_Degree_Input = 5




!   Set the Limit on the Spherical Harmonic Expansion !
L_Limit_Input = 0




!   Set the Inner Radius of the Problem !
Inner_Radius = 0.0_idp




!   Set the Outer Radius of the Problem !
Outer_Radius = 1.0_idp




!   Set the Number of Radial Elements   !
R_Elements_Input = 100




!   Set the Number of Theta Elements    !
T_Elements_Input = 10                   !   For a 1-Dimensional Simulation Set equal to 1




!   Set the Number of Phi Elements      !
P_Elements_Input = 10                   !   For a 1 or 2-Dimensional Simulations Set equal to 1





!
!   A user defines a mesh in Poseidon by providing a vector per active dimension    !
!   with each value in the vector giving the length of the associated element.      !
!   The sum of the values in the R vector should equal R_Outer - R_Inner, while     !
!   the theta and phi vectors should sum to pi and 2*pi respectively.               !
!                                                                                   !
ALLOCATE(Input_Delta_R_Vector(1:R_Elements_Input))
ALLOCATE(Input_Delta_T_Vector(1:T_Elements_Input))
ALLOCATE(Input_Delta_P_Vector(1:P_Elements_Input))


!                                               !
!   For this example, we build uniform meshes.  !
!                                               !
Input_Delta_R_Vector = (Outer_Radius - Inner_Radius)/REAL(R_Elements_Input)
Input_Delta_T_Vector = (pi)/REAL(T_Elements_Input)
Input_Delta_P_Vector = (2*pi)/REAL(P_Elements_Input)









!                           !
!   Source Input Variables  !
!                           !

!
!   For a Newtonian simulation, the density must be provided as source input.   !
!   To input to Poseidon, the user provides density value at locations within   !
!   each element.  To relate this information, the user specifies a x-space by  !
!   providing left and right limits.  Then the user specifices how many input   !
!   nodes per dimension will be provided.  For reduced dimension simulations    !
!   set the number of theta and phi nodes to 1. Next the user specifies the     !
!   quadrature points in the x-space defined by the limits for each dimension.  !
!   If a dimension is being ignored then the quadrature does not need to be     !
!   specified.                                                                  !
!                                                                               !
Left_Limit = -0.5_idp
Right_Limit = 0.5_idp

Num_Input_Nodes(1) = 5
Num_Input_Nodes(2) = 1
Num_Input_Nodes(3) = 1


!
!   The number of degrees of freedom (DOF) per element is the product   !
!   of the number of nodes per dimension.   This value represents the   !
!   total number of nodes per element where the input will be specified.!
!   This value is needed to allocate space for input variable, RHO      !
!                                                                       !
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
































!                       !
!   Output Variables    !
!                       !


!
!   I provide two kinds of output methods at the moment.  The first !
!   provides the potential at a single specific location.  The      !
!   other provides the potential in a per element fashion like the  !
!   source input above.                                             !
!                                                                   !

!                                                                           !
!   This value is used in a loop that demonstrates single location output   !
!                                                                           !
NUM_SAMPLES = 20




!
!   For multiple location output, the subroutine "Poseidon_Newtonian_Potential_Output" is called    !
!   This subroutine take input of the exact same kind as "Poseidon_Newtonian_Source_Input".  In     !
!   fact, a user can reuse the exact same input to have Poseidon return the solution values at the  !
!   source value locations.  On the other hand, a user could specify a complete new set of          !
!   locations to recieve the solution.                                                              !
!                                                                                                   !

Output_Left_Limit = -0.5
Output_Right_Limit = 0.5

Num_Output_Nodes(1) =  3
Num_Output_Nodes(2) =  1
Num_Output_Nodes(3) =  1

Num_Out_DOF = Num_Output_Nodes(1)*Num_Output_Nodes(2)*Num_Output_Nodes(3)

ALLOCATE(Potential_Output(1:Num_Out_DOF, 1:R_Elements_Input, 1:T_Elements_Input, 1:P_Elements_Input))
ALLOCATE(Output_R_Quad(1:Num_Output_Nodes(1)), Output_T_Quad(1:Num_Output_Nodes(2)), Output_P_Quad(1:Num_Output_Nodes(3)))



!                                                                           !
!   We want to recieve the solution at the Lobatto locations again, but     !
!   these may not be the same Lobatto points as before.                     !
!                                                                           !
Output_R_Quad = Initialize_LG_Quadrature_Locations(Num_Output_Nodes(1))
Output_T_Quad = Initialize_LG_Quadrature_Locations(Num_Output_Nodes(2))
Output_P_Quad = Initialize_LG_Quadrature_Locations(Num_Output_Nodes(3))




!                                                                           !
!   As with the Source input, we need to map the Lobatto quadrature points  !
!   input the space bracketed by Output_Left_Limit and Output_Right_Limit   !
!                                                                           !
Output_R_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Output_R_Quad)
Output_T_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Output_T_Quad)
Output_P_Quad = Map_From_X_Space(Left_Limit, Right_Limit, Output_P_Quad)















!                                         !
!!                                       !!
!!!         Initialize Poseidon         !!!
!!                                       !!
!                                         !


!
!   Here we truely begin to run Poseidon.  To begin it must be initialized.  There are  !
!   two methods that this can be done. 
!
!   The first method uses a single call that sets key variables for Poseidon as well    !
!       as defines the spacial meshes that will be used.  This is convienient if the    !
!       mesh will remain static throughout the lifetime of the simulation.              !
!
!   The second method uses two calls.  The first sets the key variables for Poseidon,   !
!       while the second call defines the spatial meshes.  This form is useful if the   !
!       user wishes to initialize Poseidon outside of a loop and define the mesh within.!
!                                                                                       !


IF (    1 == 1  ) THEN  !Single Call Method

    CALL Poseidon_Initialize(   FEM_Degree_Input,       &
                                L_Limit_Input,          &
                                Inner_Radius,           &
                                Outer_Radius,           &
                                R_Elements_Input,       &
                                T_Elements_Input,       &
                                P_Elements_Input,       &
                                Input_Delta_R_Vector)!, &       ! The following two input vectors are optional.
!                               Input_Delta_T_Vector,   &       ! Does not need to be included for 1D problems
!                               Input_Delta_P_Vector         )  ! Does not need to be includes for 1D or 2D problems


ELSE    ! Two Call Method

    CALL Poseidon_Initialize(   FEM_Degree_Input,       &
                                L_Limit_Input,          &
                                Inner_Radius,           &
                                Outer_Radius,           &
                                R_Elements_Input,       &
                                R_Elements_Input,       &
                                P_Elements_Input           )


    CALL Poseidon_Set_Mesh(     Input_Delta_R_Vector)!,    &     ! The following two input vectors are optional.
!                                Input_Delta_T_Vector,     &       ! Does not need to be included for 1D problems
!                                Input_Delta_P_Vector         )    ! Does not need to be includes for 1D or 2D problem




END IF





















!                                         !
!!                                       !!
!!!         Create Source Values        !!!
!!                                       !!
!                                         !


!                                                                       !
!   This part will normally be done with the user's code, but for this  !
!       example we create input data with the "Initialize_Test_Probelm" !
!       subroutine.  This routine uses a number to select which         !
!       problem is being run.  Then it creates the right source data,   !
!       Rho, and associates the pointer function, Analytic Solution,    !
!       to the right function containing the solution to the selected   !
!       problem.                                                        !
!                                                                       !





    !!!  Call To Run Test Problems !!!
CALL Poseidon_Initialize_Test_Problem(1, Num_Input_Nodes, Rho)


!                               !
!   Boundary Value Variables    !
!                               !

!
!   As it is currently, Test type 1 has this as an enclosed mass.  This can be used to
!   specify a Neumann Boundary condition.
!

Enclosed_Mass = (4.0_idp/3.0_idp)*pi!*0.5*0.5*0.5



!                                         !
!!                                       !!
!!!       Set Boundary Conditions       !!!
!!                                       !!
!                                         !


!                                                                           !
!   For the Poisson equation to be solved uniquely two boundary conditions  !
!   must be specified.  Currently Poseidon can handle Dirichlet and Neumann !
!   boundary conditions, with the hope of adding Robin boundary conditions  !
!   soon.  To set the boundary conditions the subroutine "Poseidon_Set_-    !
!   Uniform_Boundary_Condition must be called. A non-uniform input is in    !
!   development.                                                            !
!                                                                           !


!
!   Poseidon_Set_Uniform_Boundary_Condtion(BC_Location_Input, BC_Type_Input, BC_Value_Input)
!
!
!   Input Variables     :
!
!           BC_Location_Input       -   "I" for Inner Boundary
!                                       "O" for Outer Boundary
!
!           BC_Type_Input           -   "D" for Dirichlet Boundary Condition
!                                       "N" for Neumann Boundary Condition
!
!           BC_Value_Input          -   For a Dirichlet Boundary Condition specify the
!                                           Newtonian potential at the boundary.
!                                       For a Neumann Boundary Condition specify the
!                                           radial derivative of the potential at the
!                                           boundary, this is equal to the enclosed mass.
!

!CALL Poseidon_Set_Uniform_Boundary_Condition("I","D",Analytic_Solution(R_INNER, 0.0_idp, 0.0_idp))
CALL Poseidon_Set_Uniform_Boundary_Condition("I","N",0.0_idp)


PRINT*," "
PRINT*,Analytic_Solution(R_OUTER, 0.0_idp, 0.0_idp), -Enclosed_Mass/R_OUTER, -1.5*Enclosed_Mass
PRINT*," "


!CALL Poseidon_Set_Uniform_Boundary_Condition("O","D",Analytic_Solution(R_OUTER, 0.0_idp, 0.0_idp))
CALL Poseidon_Set_Uniform_Boundary_Condition("O","D",-Enclosed_Mass/R_OUTER)




!CALL Poseidon_Set_Uniform_Boundary_Condition("O","N",Enclosed_Mass/(R_OUTER*R_OUTER))













!                                         !
!!                                       !!
!!!        Input Source Values          !!!
!!                                       !!
!                                         !



!
!   This is where the source data is set.  This will need to be called everytime the source     !
!   changes. The subroutine "Poseidon_Newtonian_Source_Input" is an overloaded subroutine       !
!   whose different forms allow for simplified input for 1, 2, and 3 dimensional sources.       !
!   The three different input forms are very similar and are exemplified below.  The function   !
!   of this subroutine is to take in data from the user and map it to Poseidon's work mesh.     !
!
!   1 - Dimensional Call
!
!CALL Poseidon_Newtonian_Source_Input(       Left_Limit,             &
!                                            Right_Limit,            &
!                                            Num_Input_Nodes(1),     &
!                                            Input_R_Quad,           &
!                                            Rho                         )
!
!   2 - Dimensional Call
!
!CALL Poseidon_Newtonian_Source_Input(       Left_Limit,             &
!                                            Right_Limit,            &
!                                            Num_Input_Nodes(1:2),   &
!                                            Input_R_Quad,           &
!                                            Input_T_Quad,           &
!                                            Rho                         )
!
!   3 - Dimensional Call
!
!
!CALL Poseidon_Newtonian_Source_Input(       Left_Limit,             &
!                                            Right_Limit,            &
!                                            Num_Input_Nodes(1:3),   &
!                                            Input_R_Quad,           &
!                                            Input_T_Quad,           &
!                                            Input_P_Quad,           &
!                                            Rho                         )
!
!
!
!   Notice difference in the size on the Num_Input_Nodes input vector, as well as
!   the additional quadrature input vectors for higher dimensions.  
!
!
!   Input Variables :
!
!   Left, Right_Limit   -   Real Values, Left and Right limits on the space in which the source
!                           locations are given.
!
!   Num_Input_Nodes           -   Integer Valued Vector, Length(1:1-3, dependent on input form -> dimensionality of input)
!                           Each value gives the number of source locations per dimension.
!
!   Input_*_Quad        -   Real Valued Vector, Each value gives the location of source data in the quadrature 
!                           space between Left_Limit and Right_Limit for the dimension *, where * = R (radial), T (theta),
!                           or P (phi). For 1-Dimensional input only Input_R_Quad need be provided, for 2-Dimensional 
!                           input only Input_R_Quad, and Input_T_Quad.
!
!   Rho                 -   Real Valued Vector, This is the main source vector which contains the density value at each
!                           point described by the Input_*_Quad vectors.




CALL Poseidon_Newtonian_Source_Input(       Left_Limit,             &
                                            Right_Limit,            &
                                            Num_Input_Nodes(1),           &
                                            Input_R_Quad,           &
                                            Rho                         )



















!                                         !
!!                                       !!
!!!             Run Poseidon            !!!
!!                                       !!
!                                         !


!                                                                                                       !
!   This subroutine actually solves the gravitational Poisson equation given the source                 !
!   provided by the user using a mixed method consisting of a spectral expansion in the                 !
!   angular dimensions using the Spherical Harmonic functions and the Finite Element method             !
!   to solve for the expansion coefficients. No further input is required at this time.                 !
!                                                                                                       !
!   Also it should be noted that currently there are limited safety-checks to make sure that the        !
!   user has provided all the necessary input for Poseidon to function properly.  This feature is       !
!   currently in development.                                                                           !
!                                                                                                       !

CALL Poseidon_Run()
















!                                                                                               !
!   There are two main methods for output at the moment. The first provides the solution at     !
!   at specific location within the computational domain.  The second provides vector output    !
!   with the values giving the solution at locations specified in the same manner as the        !
!   source is input into Poseidon.  Below we give demonstrations of both.                       !
!                                                                                               !

!
!   Specific Locaiton Output
!
! This can be toggled on and off here.
IF ( 1 == 1 ) THEN



    ! This helps us build output locations.
    deltar = (R_OUTER - R_INNER)/REAL(NUM_SAMPLES, KIND = idp)




    DO i = 0,NUM_SAMPLES

        r = i*deltar + R_INNER
        theta = pi/2.0_idp
        phi = pi/2.0_idp




        !
        !   As with many Poseidon calls, there are multiple versions for different dimensions.
        !   The function call differs here only in the number of dimensions given to specify the 
        !   location.  If a value is not provided it is chosen to be 0.0_idp.
        !


        !CALL potential = Calculate_Potential_At_Location( r, theta, phi )
        !CALL potential = Calculate_Potential_At_Location( r, theta )
        potential =  Calculate_Potential_At_Location( r )












        !   By calling "Poseidon_Initialize_Test_Problem" we have acess to an analytic solution
        !   for the chosen problem through the function Analytic_Solution. This function currently
        !   requires all three values of input even for 1 and 2 dimensional problems.  This will be 
        !   fixed in a future update.
        TMP = Analytic_Solution(r,theta,phi)




        !
        !   This compares Poseidon's solution with the analytic solution and prints the results.
        !   The if statement allows for exact solutions by avoiding a division by zero.
        !

        IF ( ABS(TMP) == 0.0_idp) THEN

            print*,r,potential,TMP, ABS(TMP-potential)/ABS(1e-16)

        ELSE

            print*,r,potential,TMP, ABS(TMP-potential)/ABS(TMP)

        END IF

    END DO

    PRINT*," "

END IF









!
!   Multiple Locaiton Output
!
! This can be toggled on and off here.
IF (1 == 0) THEN




!   For multiple location output in a single call, the subroutine "Poseidon_Newtonian_Potential_Output"
!   can be used.  This subroutine functions much like "Poseidon_Newtonian_Source_Input" with the same 
!   call forms for different dimensions.  The difference between the two being that 
!   "Poseidon_Newtonian_Potential_Output" takes an empty vector and fills it with solution values for 
!   the locations described by the other input variables.
!
!   1 - Dimensional Call
!
!CALL Poseidon_Newtonian_Potential_Output(      Left_Limit,             &
!                                               Right_Limit,            &
!                                               Num_Output_Nodes(1),    &
!                                               Input_R_Quad,           &
!                                               Potential                   )
!
!   2 - Dimensional Call
!
!CALL Poseidon_Newtonian_Potential_Output(      Left_Limit,             &
!                                               Right_Limit,            &
!                                               Num_Output_Nodes(1:2),  &
!                                               Output_R_Quad,          &
!                                               Output_T_Quad,          &
!                                               Potential                   )
!
!   3 - Dimensional Call
!
!
!CALL Poseidon_Newtonian_Potential_Output(      Left_Limit,             &
!                                               Right_Limit,            &
!                                               Num_Output_Nodes(1:3),  &
!                                               Output_R_Quad,          &
!                                               Output_T_Quad,          &
!                                               Output_P_Quad,          &
!                                               Potential                   )
!
!
!
!   Notice difference in the size on the Num_Input_Nodes input vector, as well as
!   the additional quadrature input vectors for higher dimensions.
!
!
!   Input Variables :
!
!   Left, Right_Limit   -   Real Values, Left and Right limits on the space in which the solution
!                           locations are desired.
!
!   Num_Input_Nodes           -   Integer Valued Vector, Length(1:1-3, dependent on output form -> dimensionality of output)
!                           Each value gives the number of solution locations per dimension.
!
!   Output_*_Quad       -   Real Valued Vector, Each value gives the location where the solution is desired in the quadrature
!                           space between Left_Limit and Right_Limit for the dimension *, where * = R (radial), T (theta),
!                           or P (phi). For 1-Dimensional input only Output_R_Quad need be provided, for 2-Dimensional
!                           input only Output_R_Quad, and Output_T_Quad.
!
!   Potential           -   Real Valued Vector, This is the vector where the solution value for each location
!                           described by the Output_*_Quad vectors will be stored at the completion of the subroutine.



    CALL Poseidon_Newtonian_Potential_Output(   Left_Limit,             &
                                                Right_Limit,            &
                                                Num_Output_Nodes(1),    &
                                                Output_R_Quad,          &
                                                Potential_Output            )


    DO i = 1,NUM_R_ELEMENTS

        PRINT*,Potential_Output(:,i,1,1)

    END DO


END IF








!                                         !
!!                                       !!
!!!             Close Poseidon          !!!
!!                                       !!
!                                         !


!
!   This subroutine closes down Poseidon by deallocating all space that was
!   allocated for Poseidon's use.
!

CALL Poseidon_Close()




















END PROGRAM Poseidon_Newtonian
