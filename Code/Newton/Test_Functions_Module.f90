   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Test_Functions_Module                                                        !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!   Contains the functions and subroutines used in the test codes provided.      !##!
!##!    These functions initialize specific sources that produce known potentials   !##!
!##!    which allow for comparison of the code's solution to the problem and the    !##!
!##!    known analytic solution.                                                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Poseidon_Initialize_Test_Problem                                    !##!
!##!                                                                                !##!
!##!    +201+   Test_Source_Spherical_Symmetry_No_Surface                           !##!
!##!    +202+   Test_Source_Spherical_Symmetry_With_Surface                         !##!
!##!    +203+   Source_MacLaurin_Ellipsoid                                          !##!
!##!                                                                                !##!
!##!    +301+   Test_Spherical_Symmetry_No_Surface                                  !##!
!##!    +302+   Test_Spherical_Symmetry_With_Surface                                !##!
!##!    +303+   Test_MacLaurin_Ellipsoid                                            !##!
!##!                                                                                !##!
!##!    +401+   MacLaurin_Radius                                                    !##!
!##!    +402+   MacLaurin_Root_Finder                                               !##!
!##!    +403+   VECTOR_PERTURB                                                      !##!
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
                ONLY :  Source_Function_Flag, POWER_A, RHO_O,               &
                        R_MIN, R_MAX, SPHEROID_TYPE, R_INNER, R_OUTER,      &
                        NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS,     &
                        INNER_BC_TYPE, OUTER_BC_TYPE,                       &
                        INNER_DIR_BC_INPUT, INNER_NEU_BC_INPUT,                         &
                        OUTER_DIR_BC_INPUT, OUTER_NEU_BC_INPUT,                         &
                        Star_Surface,                                       &
                        Source_Function, Source_Degrees,                    &
                        Source_Term_Coefficients,                           &
                        Analytic_Solution,                                  &
                        rlocs, tlocs, plocs,                                &
                        Test_Source_Input, Test_Space_Allocated_Flag,       &
                        Test_Run_Flag


USE Additional_Functions_Module, &
                ONLY :  Lagrange_Poly, Map_From_X_Space,                   &
                        Initialize_LG_Quadrature,                           &
                        Initialize_LGL_Quadrature_Locations,                &
                        Initialize_LG_Quadrature_Locations,                 &
                        MVMULT_FULL







IMPLICIT NONE











!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS

 !+101+############################################################################!
!                                                                                   !
!       Poseidon_Initialize_Test_Problem                                            !
!                                                                                   !
!      This code initializes the source for the selected problem as well as points  !
!       the code to the analytic solution associated with the given source.         !
!                                                                                   !
!-----------------------------------------------------------------------------------!
!                                                                                   !
!   Input                                                                           !
!                                                                                   !
!           Test_Number     -   Integer, Selects problem to be solved.              !
!                                                                                   !
!                                   1)  Spherically Symmetric Source with radial    !
!                                           dependence.                             !
!                                                                                   !
!                                   2) Spherically Symmetric source with radial     !
!                                       dependence and include surface discontinuity!
!                                                                                   !
!                                   3) MacLaurin Ellipsoid                          !
!                                                                                   !
 !#################################################################################!
SUBROUTINE Poseidon_Initialize_Test_Problem(Test_Number, Num_Nodes, Test_Source_Input, Input_Delta_R_Vector)

INTEGER,                        INTENT(IN)                                  ::  Test_Number
INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes



REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:NUM_R_ELEMENTS,                           &
                            1:NUM_T_ELEMENTS,                           &
                            1:NUM_P_ELEMENTS),                          &
                            INTENT(INOUT)                                   ::  Test_Source_Input



REAL(KIND = idp), OPTIONAL, DIMENSION(1:NUM_R_ELEMENTS), INTENT(INOUT)      :: Input_Delta_R_Vector





INTEGER                                                                     ::  Num_DOF



Test_Run_Flag = .TRUE.


Num_DOF = Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3)



INNER_BC_TYPE = "DRCH"
OUTER_BC_TYPE = "DRCH"




IF (Test_Number .EQ. 1) THEN


    Source_Function_Flag = 1


    STAR_SURFACE = R_OUTER + 1.0_idp



    CALL Test_Source_Spherical_Symmetry(Num_Nodes, Test_Source_Input)

    Analytic_Solution => Test_Spherical_Symmetry_No_Surface







ELSE IF (Test_Number .EQ. 2) THEN


    Source_Function_Flag = 2


    Star_Surface = 0.50_idp*(R_OUTER - R_INNER)+ R_INNER


    CALL Test_Source_Spherical_Symmetry(Num_Nodes, Test_Source_Input)

    POWER_A = 0
    Analytic_Solution => Test_Spherical_Symmetry_With_Surface






ELSE IF (Test_Number .EQ. 3) THEN

    Source_Function_Flag = 3

    STAR_SURFACE = R_OUTER + 1.0_idp

    R_MIN = 0.10_idp

    R_MAX = 0.25_idp




    Analytic_Solution => Test_MacLaurin_Ellipsoid







    !!!     To aid the integration of the source we want design a mesh were the     !!!
    !!!     surface discontinuity of the star is contained within a single cell.    !!!
    !!!     This creates a Delta_R_Vector that can be used in Poseidon that         !!!
    !!!     creates such a mesh.                                                    !!!

    !!!     1 Radial Element        !!!
    IF ( NUM_R_ELEMENTS == 1 ) THEN



        Input_Delta_R_Vector(1) = R_OUTER - R_INNER



    !****!     Odd Number of Radial Elements   !****!
    ELSE IF ( MOD(NUM_R_ELEMENTS,2) == 1) THEN






        !!!     Fill Width of Elements Before Discontinuity     !!!
        Input_Delta_R_Vector( 1:(NUM_R_ELEMENTS-1)/2 ) = (R_MIN - R_INNER)/((NUM_R_ELEMENTS - 1)/2)


        !!!     Central element entirely contains the star surface  !!!
        Input_Delta_R_Vector( (NUM_R_ELEMENTS + 1)/2 ) = R_MAX - R_MIN



        !!!     Fill Width of Elements After Discontinuity      !!!
        Input_Delta_R_Vector(((NUM_R_ELEMENTS + 1)/2)+1:NUM_R_ELEMENTS) = (R_OUTER - R_MAX)/((NUM_R_ELEMENTS - 1)/2)








    !****!     Even Number of Radial Elements   !****!
    ELSE


        !!!     Fill Width of Elements Before Discontinuity     !!!
        Input_Delta_R_Vector( 1:(NUM_R_ELEMENTS/2)-1 ) = (R_MIN - R_INNER)/((NUM_R_ELEMENTS/2)-1)


        !!! Find the width of elements not containing the discontinuity !!!
        Input_Delta_R_Vector( NUM_R_ELEMENTS/2 ) = R_MAX - R_MIN



        !!!     Fill Width of Elements After Discontinuity      !!!
        Input_Delta_R_Vector((NUM_R_ELEMENTS/2)+1:NUM_R_ELEMENTS) = (R_OUTER - R_MAX)/(NUM_R_ELEMENTS/2)




    END IF







END IF


END SUBROUTINE Poseidon_Initialize_Test_Problem

















!+201+######################
!
!   Spherically Symmetric
!
!###########################
SUBROUTINE Test_Source_Spherical_Symmetry(Num_Nodes, Test_Source_Input)



INTEGER,    DIMENSION(1:3),     INTENT(IN)                                  ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:NUM_R_ELEMENTS,                           &
                            1:NUM_T_ELEMENTS,                           &
                            1:NUM_P_ELEMENTS),                          &
                            INTENT(INOUT)                                   ::  Test_Source_Input

                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !
INTEGER                                                     ::  re, te, pe, i,                  &
                                                                Output_R, Output_T, Output_P,   &
                                                                l, m, d,                        &
                                                                Output_Here,                    &
                                                                Num_Output_DOF


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_X_Locations


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_Locations

REAL(KIND = idp)                                            ::  TMP





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
                        !!!     Map Output Locations to [ -1, 1 ] Space      !!!
                         !!                                                 !!
                          !                                                 !

Output_R_X_Locations = Initialize_LG_Quadrature_Locations(Num_Nodes(1))






DO re = 0,NUM_R_ELEMENTS-1


     !                                                          !
    !!   Map Radial Locations from [-1,1] space to real space.   !!
     !                                                          !
    Output_R_Locations = Map_From_X_Space(rlocs(re), rlocs(re + 1), Output_R_X_Locations)



    DO te = 0,NUM_T_ELEMENTS-1



        DO pe = 0,NUM_P_ELEMENTS-1


             !                                          !
            !!   Set/Reset Output Vector Location       !!
             !                                          !
            Output_Here = 1





            DO Output_P = 1,Num_Nodes(3)

                DO Output_T = 1,Num_Nodes(2)

                    DO Output_R = 1,Num_Nodes(1)


                        TMP = 1.0_idp







                        IF (Output_R_Locations(Output_R) .LE. Star_Surface) THEN


                            IF (POWER_A >= 0) THEN

                                DO i = 1,POWER_A

                                    TMP = TMP*Output_R_Locations(Output_R)

                                END DO
                            END IF


                            IF (POWER_A < 0) THEN
                                DO i = 1,ABS(POWER_A)

                                    TMP = TMP/Output_R_Locations(Output_R)

                                END DO
                            END IF


                            Test_Source_Input(Output_Here, re+1, te+1, pe+1) = RHO_O*TMP


                        ELSE

                            Test_Source_Input(Output_Here, re+1, te+1, pe+1) = 0.0_idp


                        END IF



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





END SUBROUTINE Test_Source_Spherical_Symmetry














!+201+######################
!
!   Test Source MacLaurin
!
!###########################
SUBROUTINE Test_Source_MacLaurin(Num_Nodes, Test_Source_Input)



INTEGER,    DIMENSION(1:3),     INTENT(IN)                  ::  Num_Nodes


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1)*Num_Nodes(2)*Num_Nodes(3),   &
                            1:NUM_R_ELEMENTS,                           &
                            1:NUM_T_ELEMENTS,                           &
                            1:NUM_P_ELEMENTS),                          &
                            INTENT(INOUT)                                   ::  Test_Source_Input

                         !                        !
                        !!  Subroutine Variables  !!
                         !                        !
INTEGER                                                     ::  re, te, pe, i,                  &
                                                                Output_R, Output_T, Output_P,   &
                                                                l, m, d,                        &
                                                                Output_Here,                    &
                                                                Num_Output_DOF


REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_X_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2))                 ::  Output_T_X_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3))                 ::  Output_P_X_Locations

REAL(KIND = idp), DIMENSION(1:Num_Nodes(1))                 ::  Output_R_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(2))                 ::  Output_T_Locations
REAL(KIND = idp), DIMENSION(1:Num_Nodes(3))                 ::  Output_P_Locations

REAL(KIND = idp)                                            ::  TMP


REAL(KIND = idp)                                            ::  A, B, C, VAL
REAL(KIND = idp)                                            ::  rsqr, sinsqr



IF (SPHEROID_TYPE .EQ. 'OBLATE') THEN

    A = R_MAX
    B = R_MAX
    C = R_MIN

ELSE IF (SPHEROID_TYPE .EQ. 'PROLATE') THEN

    A = R_MAX
    B = R_MIN
    C = R_MIN

END IF


A = A*A
B = B*B
C = C*C




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
                        !!!     Map Output Locations to [ -1, 1 ] Space      !!!
                         !!                                                 !!
                          !                                                 !

Output_R_X_Locations = Initialize_LG_Quadrature_Locations(Num_Nodes(1))
Output_T_X_Locations = Initialize_LG_Quadrature_Locations(Num_Nodes(2))
Output_P_X_Locations = Initialize_LG_Quadrature_Locations(Num_Nodes(3))






DO re = 0,NUM_R_ELEMENTS-1


     !                                                          !
    !!   Map Radial Locations from [-1,1] space to real space.  !!
     !                                                          !
    Output_R_Locations = Map_From_X_Space(rlocs(re), rlocs(re + 1), Output_R_X_Locations)



    DO te = 0,NUM_T_ELEMENTS-1


         !                                                          !
        !!   Map Theta Locations from [-1,1] space to real space.   !!
         !                                                          !
        Output_T_Locations = Map_From_X_Space(tlocs(te), tlocs(te + 1), Output_T_X_Locations)


        DO pe = 0,NUM_P_ELEMENTS-1


             !                                                          !
            !!   Map Phi Locations from [-1,1] space to real space.  !!
             !                                                          !
            Output_P_Locations = Map_From_X_Space(plocs(pe), plocs(pe + 1), Output_P_X_Locations)

             !                                          !
            !!   Set/Reset Output Vector Location       !!
             !                                          !
            Output_Here = 1





            DO Output_P = 1,Num_Nodes(3)

                DO Output_T = 1,Num_Nodes(2)

                    DO Output_R = 1,Num_Nodes(1)




                        rsqr = Output_R_Locations(Output_R)*Output_R_Locations(Output_R)
                        sinsqr = sin(Output_T_Locations(Output_T))*sin(Output_T_Locations(Output_T))

                        VAL = (rsqr*cos(Output_P_Locations(Output_P))*cos(Output_P_Locations(Output_P))*sinsqr)/A   &
                            + (rsqr*sin(Output_P_Locations(Output_P))*sin(Output_P_Locations(Output_P))*sinsqr)/B   &
                            + (rsqr*cos(Output_T_Locations(Output_T))*cos(Output_T_Locations(Output_T)))/C


                        IF (VAL .LE. 1.0_idp) THEN

                            Test_Source_Input(Output_Here, re, te, pe) = RHO_O

                        ELSE

                            Test_Source_Input(Output_Here, re, te, pe) = 0.0_idp

                        END IF





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





END SUBROUTINE Test_Source_MacLaurin




















!+301+##################################################################
!
!   Test_Spherical_Symmetry_No_Surface
!
!#######################################################################
PURE FUNCTION Test_Spherical_Symmetry_No_Surface(r, theta, phi)

REAL(KIND = idp),INTENT(IN)         :: r, theta, phi
REAL(KIND = idp)                    :: Test_Spherical_Symmetry_No_Surface

INTEGER                             :: i
REAL(KIND = idp)                    :: TMP





TMP = 1.0_idp
IF (POWER_A+2 > 0) THEN
    DO i = 1,POWER_A+2

        TMP = TMP*r

    END DO



    Test_Spherical_Symmetry_No_Surface = -(4.0_idp/3.0_idp)* pi *RHO_O * (3*R_OUTER*R_OUTER - r*r)/(2*R_OUTER*R_OUTER*R_OUTER)

    !Test_Spherical_Symmetry_No_Surface =  4.0_idp * pi *RHO_O*TMP/((3+POWER_A)*(2+POWER_A))

END IF






IF (POWER_A + 2 .EQ. 0) THEN

    Test_Spherical_Symmetry_No_Surface = -4.0_idp * pi *RHO_O * log(r)

END IF







IF (POWER_A + 2 .EQ. -1) THEN

    Test_Spherical_Symmetry_No_Surface = 4.0_idp * pi *RHO_O*(log(r)+1)/r

END IF








IF (POWER_A+2 < -1) THEN
    DO i = 1,ABS(POWER_A + 2)

        TMP = TMP/r

    END DO

    Test_Spherical_Symmetry_No_Surface = -4.0_idp * pi *RHO_O*TMP/((3+POWER_A)*(2+POWER_A))

END IF





END FUNCTION Test_Spherical_Symmetry_No_Surface













!+302+#################################################################
!
!   Test_Spherical_Symmetry_With_Surface
!
!#######################################################################
FUNCTION Test_Spherical_Symmetry_With_Surface(r, theta, phi)

REAL(KIND = idp),INTENT(IN)         :: r, theta, phi
REAL(KIND = idp)                    :: Test_Spherical_Symmetry_With_Surface

INTEGER                             :: i
REAL(KIND = idp)                    :: TMP_r, TMP_BIG_R




IF (r .LE. STAR_SURFACE) THEN


    TMP_r = 1.0_idp
    TMP_BIG_R = 1.0_idp





    IF (POWER_A+2 > 0) THEN
        DO i = 1,POWER_A+2

            TMP_r = TMP_r*r
            TMP_BIG_R = TMP_BIG_R*STAR_SURFACE


        END DO

        Test_Spherical_Symmetry_With_Surface = (1.0_idp/((POWER_A+3.0_idp)*(POWER_A + 2.0_idp)))   &
                                                * 4.0_idp * pi *RHO_O                               &
                                                * (TMP_r - (POWER_A + 3)*TMP_BIG_R)


    END IF






    IF (POWER_A + 2 .EQ. 0) THEN

        Test_Spherical_Symmetry_With_Surface = 4.0_idp * pi *RHO_O * log(r)

    END IF







    IF (POWER_A + 2 .EQ. -1) THEN

        Test_Spherical_Symmetry_With_Surface = 4.0_idp * pi *RHO_O*(log(r)+1)/r

    END IF








    IF (POWER_A+2 < -1) THEN
        DO i = 1,ABS(POWER_A + 2)

            TMP_r = TMP_r/r

        END DO

        Test_Spherical_Symmetry_With_Surface = RHO_O*TMP_r/((3+POWER_A)*(2+POWER_A))

    END IF


ELSE  ! r > STAR_SURFACE


    Test_Spherical_Symmetry_With_Surface = -(4.0_idp/3.0_idp) * pi *RHO_O                           &
                                            *(STAR_SURFACE*STAR_SURFACE*STAR_SURFACE)/r

END IF







END FUNCTION Test_Spherical_Symmetry_With_Surface










!+303+##########################################################
!
!   MacLaurin_Potential - Calulates the potential at a given location for the MacLaurin Spheriod
!
!################################################################
PURE FUNCTION Test_MacLaurin_Ellipsoid (r, theta, phi)


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi


REAL(KIND = idp)                                            ::  Test_MacLaurin_Ellipsoid

REAL(KIND = idp)                                            ::  A, B, C, X, Y, Z, VAL, &
                                                                rsqr, psqr, &
                                                                ecc, eccsqr, eccmod, &
                                                                Parta, Partb, Partc, &
                                                                lambda, POTENTIAL,h


!
! Determine if the location is inside or outside of the spheroid
!
IF (SPHEROID_TYPE .EQ. 'OBLATE') THEN

    A = R_MAX
    B = R_MAX
    C = R_MIN

ELSE IF (SPHEROID_TYPE .EQ. 'PROLATE') THEN

    A = R_MAX
    B = R_MIN
    C = R_MIN

END IF



!   Calculate squares of semi-axes
A = A*A
B = B*B
C = C*C



!   Calculate Eccentricity !
ecc =   sqrt(1 - C/A)
eccsqr = ecc*ecc




!   Calculate Cartesian Coordinates
X = r * cos(phi) * sin(theta)
Y = r * sin(phi) * sin(theta)
Z = r * cos(theta)



VAL = (X*X)/A + (Y*Y)/B + (Z*Z)/C
!
!   VAL > 1 -> point is outside.  VAL <= 1 -> point is inside
!
IF (VAL .LE. 1.0_idp) THEN

    !
    !   Equations for potential inside homegenous spheroid
    !
    !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3
    !
    IF (SPHEROID_TYPE .EQ. 'OBLATE') THEN

        eccmod = sqrt(1- eccsqr)/(eccsqr*ecc)*ASIN(ecc)

        Parta = eccmod - (1 - eccsqr)/eccsqr

        Partb = 2/eccsqr - 2*eccmod



        POTENTIAL = pi*((2*A - X*X - Y*Y)*Parta + (C - Z*Z)*Partb)


    ELSE IF (SPHEROID_TYPE .EQ. 'PROLATE') THEN

        eccmod = (1 - eccsqr)/(eccsqr*ecc) * LOG((1+ecc)/(1-ecc))

        Parta = eccmod - 2*(1 - eccsqr)/eccsqr

        Partb = 1/eccsqr - 0.5_idp*eccmod


        POTENTIAL = pi*(2*C - Y*Y - Z*Z)*Partb + (A - X*X)*Parta


    END IF


ELSE IF (VAL > 1.0_idp) THEN


    !
    !   Equations for potential outside homegenous spheroid
    !
    !       REF:    Ellipsoidal Figures of Equilibrium - Chandrasekhar - Chapter 3 -
    !
    lambda = MacLaurin_Root_Finder(r, theta, phi)

    Parta = ATAN(sqrt((R_MAX*R_MAX - R_MIN*R_MIN)/(R_MIN*R_MIN + lambda)))

    Partb = 1/sqrt(R_MAX*R_MAX - R_MIN*R_MIN)

    Partc = Partb*Partb*Partb





    IF (SPHEROID_TYPE .EQ. 'OBLATE') THEN




        POTENTIAL = pi*R_MAX*R_MAX*R_MIN * ( 2 * Partb * Parta &
                                    - (X*X + Y*Y)*Partc*(Parta - (sqrt(A-C)*sqrt(C + lambda))/(A + lambda)) &
                                    + 2*(Z*Z)*Partc*(Parta - sqrt(A - C)/sqrt(C + lambda)))


    ELSE IF (SPHEROID_TYPE .EQ. 'PROLATE') THEN


        POTENTIAL = R_MAX*R_MIN*R_MIN*(-Parta - Partb &
                        + (X*X)*((1.5*R_MAX - 0.5*R_MIN + lambda)*Partb*Partb + Parta/(R_MAX - R_MIN)) &
                        + (Y*Y + Z*Z) * ((1.5*R_MIN - 0.5*R_MAX + lambda)*Partc*Partc + Parta/(R_MAX - R_MIN)))



    END IF


END IF

Test_MacLaurin_Ellipsoid  = POTENTIAL


END FUNCTION Test_MacLaurin_Ellipsoid 





















!+401+######################
!
!   MacLaurin_Radius
!
!###########################
FUNCTION MacLaurin_Radius(theta, phi)


REAL(KIND = idp), INTENT(IN)            ::  theta, phi


REAL(KIND = idp)                        ::  A, B, C, RADIUS_HERE, TMP_RAD


REAL(KIND = idp)                        ::  MacLaurin_Radius



IF (SPHEROID_TYPE .EQ. 'OBLATE') THEN

    A = R_MAX
    B = R_MAX
    C = R_MIN



    RADIUS_HERE = A*C / sqrt(C*C*sin(theta)*sin(theta) + A*A*cos(theta)*cos(theta))

ELSE IF (SPHEROID_TYPE .EQ. 'PROLATE') THEN

    A = R_MAX
    B = R_MIN
    C = R_MIN



    RADIUS_HERE = A*B*SQRT(B*B*cos(theta)*cos(theta)*sin(phi)*sin(phi) &
                            + A*A*(sin(phi)*sin(phi)*sin(theta)*sin(theta) + cos(theta)*cos(theta)))


END IF




MacLaurin_Radius = RADIUS_HERE





END FUNCTION MacLaurin_Radius














!+402+#################################################################
!
!   MacLaurin_Root_Finder
!
!#######################################################################
PURE FUNCTION MacLaurin_Root_Finder(r, theta, phi)


REAL(KIND = idp), INTENT(IN)                                ::  r, theta, phi

REAL(KIND = idp)                                            ::  MacLaurin_Root_Finder

REAL(KIND = idp)                                            ::  rsqr, xsqr, ysqr, zsqr, &
                                                                asqr, csqr, B, C

!! r^2  !!
rsqr = r * r


!! x^2 = r^2 * cos(theta)^2 * sin(phi)^2 !!
xsqr = rsqr * cos(phi) * cos(phi)* sin(theta) * sin(theta)

!! y^2 = r^2 * sin(theta)^2 * sin(phi)^2 !!
ysqr = rsqr * sin(phi) * sin(phi)* sin(theta) * sin(theta)

!! z^2 = r^2 * cos(phi)^2  !!
zsqr = rsqr * cos(theta) * cos(theta)



asqr = R_MAX * R_MAX

csqr = R_MIN * R_MIN

B = asqr + csqr - xsqr - ysqr - zsqr


IF (SPHEROID_TYPE .EQ. "OBLATE") THEN


    C = asqr*(csqr - zsqr) - csqr*(xsqr + ysqr)


ELSE IF ( SPHEROID_TYPE .EQ. "PROLATE" ) THEN


    C = csqr*(asqr - xsqr) - asqr*(ysqr + zsqr)


END IF


MacLaurin_Root_Finder = 0.5 * ( - B + sqrt(B*B - 4*C))




END FUNCTION MacLaurin_Root_Finder






!+403+#################################################################!
!
!   VECTOR_PERTURB - Scales elements in a vector, VECT, randomly by
!                     a value in [-SCALE_FACTOR,+SCALE_FACTOR]
!
!######################################################################!
FUNCTION VECTOR_PERTURB(N, SCALE_FACTOR, VECT)

INTEGER, INTENT(IN)                                     ::  N
REAL(KIND = idp), INTENT(IN)                            ::  SCALE_FACTOR
!REAL(KIND = idp), DIMENSION(0:N-1), INTENT(IN)          ::  VECT
COMPLEX(KIND = idp), DIMENSION(0:N-1), INTENT(IN)          ::  VECT


REAL(KIND = idp), DIMENSION(0:N-1)                      ::  VECTOR_PERTURB


INTEGER                                                 ::  i, k
REAL(KIND = idp)                                        ::  RN, SCALER


INTEGER                                                 ::  values(1:8) !   ^
INTEGER, ALLOCATABLE, DIMENSION(:)                      ::  seed        !   |
CALL date_and_time(values = values)                                     !   |
CALL RANDOM_SEED(size = k)                                              !   |=>  Used to initalize RNG
ALLOCATE(seed(1:k))                                                     !   |
seed(:) = values(8)                                                     !   |
CALL RANDOM_SEED(put = seed)                                            !   v






DO i = 0,N-1

    CALL RANDOM_NUMBER(RN)

    RN = 2*RN - 1

    SCALER = 1 + (SCALE_FACTOR * RN * VECT(i))


    VECTOR_PERTURB(i) = SCALER * VECT(i)


END DO





END FUNCTION VECTOR_PERTURB















END MODULE Test_Functions_Module