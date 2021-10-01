   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Integrator_3D_Module                                                  !##!
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
            ONLY :  idp, pi, twopi


USE Global_Variables_And_Parameters, &
            ONLY :  NUM_R_ELEMENTS, NUM_P_ELEMENTS, NUM_T_ELEMENTS,     &
                    rlocs, plocs, tlocs,                                &
                    Source_Function_Flag, NUM_R_NODES, L_LIMIT, DEGREE, &
                    Source_Vector, Source_Function, Discontinuity_Flag, &
                    Source_Degrees, Source_Term_Coefficients

USE Additional_Functions_Module, &
            ONLY :  Legendre_Poly, Lagrange_Poly,   &
                    Norm_Factor, POWER, Map_From_X_Space, Map_To_X_Space,          &
                    Initialize_LG_Quadrature, Initialize_LGL_Quadrature

IMPLICIT NONE


CONTAINS



 !+201b+####################################################################!
!                                                                           !
!                            Triple Integral                                !
!                                                                           !
!---------------------------------------------------------------------------!
!                                                                           !
!   This function performs the triple intergral over the source function    !
!       used to defined values in the source vector of the linear system    !
!       established by the mixed spectral/fintie element method.            !
!                                                                           !
 !#########################################################################!
SUBROUTINE Triple_Integral(Src_Array)


COMPLEX(KIND = idp), DIMENSION(0:NUM_R_NODES-1, -L_LIMIT:L_LIMIT, 0:L_LIMIT), INTENT(INOUT)     ::  Src_Array


INTEGER                                             ::  l, m, re, rd, te, td, pe, pd, p


REAL(KIND = idp)                                    ::  deltar, deltap, deltat,     &
                                                        SphereHarm_NormFactor

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  R_xlocs, R_locs, R_weights, &
                                                        P_xlocs, P_locs, P_weights, &
                                                        T_xlocs, T_locs, T_weights, &
                                                        SphereHarm_ThetaPart


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_xlocs, Poly_weights, LagP


REAL(KIND = idp)                                    ::  R_SQR, T_PRE
COMPLEX(KIND = idp)                                 ::  P_PRE


  !                                         !
 !!                                         !!
!!!     Initialize Source Vector To Zero    !!!
 !!                                         !!
  !                                         !
Src_Array = 0.0_idp










  !                                                         !
 !!                                                         !!
!!!     Allocate Space for Quadratures Points and Weights   !!!
 !!                                                         !!
  !                                                         !
ALLOCATE(R_xlocs(1:Source_Degrees(1)), R_locs(1:Source_Degrees(1)), R_weights(1:Source_Degrees(1)))
ALLOCATE(P_xlocs(1:Source_Degrees(2)), P_locs(1:Source_Degrees(2)), P_weights(1:Source_Degrees(2)))
ALLOCATE(T_xlocs(1:Source_Degrees(3)), T_locs(1:Source_Degrees(3)), T_weights(1:Source_Degrees(3)))

!ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
!ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))
ALLOCATE(SphereHarm_ThetaPart(0:Source_Degrees(2)))





  !                                                 !
 !!                                                 !!
!!!     Initalize Quadrature Points and Weights     !!!
 !!                                                 !!
  !                                                 !
CALL Initialize_LGL_Quadrature(DEGREE, Poly_xlocs, Poly_weights)

CALL Initialize_LG_Quadrature(Source_Degrees(1), R_xlocs, R_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(2), P_xlocs, P_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(3), T_xlocs, T_weights)



   !                                                     !
  !!                                                     !!
 !!!     Theta and Phi integrals are done in one cell,   !!!
!!!!     so we map quadrature points from (-1,1) space   !!!!
 !!!     to theta, phi space now.                        !!!
  !!                                                     !!
   !                                                     !

P_locs = Map_From_X_Space(0.0_idp, 2*pi, P_xlocs)
T_locs = Map_From_X_Space(0.0_idp, pi, T_xlocs)

deltap = 2.0_idp * pi
deltat = pi





 !                                     !
!!                                     !!
!!     Begin Loop over all l,m pairs   !!
!!                                     !!
 !                                     !
DO l = 0, L_LIMIT

    DO m = -l, l


          !                                                                                     !
         !!                                                                                     !!
        !!!     The Spherical Harmonic functions can be subdivided into three parts:            !!!
        !!!         1)  A normialization factor dependent on l and m.                           !!!
        !!!         2)  An Associated Legendre Polynomial, dependent on l, m, and theta         !!!
        !!!         3)  An exponential part, dependent on m, and phi                            !!!
        !!!     To save repetitious calculation these three parts are calculated seperately     !!!
        !!!     and combined at a later point to recreate the whole spherical harmonic function.!!!
        !!!                                                                                     !!!
        !!!     Here we compute the first part, the normalization factor.                       !!!
         !!                                                                                     !!
          !                                                                                     !
        SphereHarm_NormFactor = POWER(-1.0_idp, m)*Norm_Factor(l,-m)




         !                              !
        !!                              !!
        !!  Loop over Radial Elements   !!
        !!                              !!
         !                              !
        DO re = 0,NUM_R_ELEMENTS - 1


              !                                                                             !
             !!                                                                             !!
            !!!     Map quadrature points from (-1,1) space to the space of element, re.    !!!
             !!                                                                             !!
              !                                                                             !
            R_locs = Map_From_X_Space(rlocs(re),rlocs(re+1), R_xlocs)
            deltar = (rlocs(re+1) - rlocs(re))



             !                              !
            !!                              !!
            !!  Loop over Theta Elements    !!
            !!                              !!
             !                              !
            DO te = 0,NUM_T_ELEMENTS - 1


                 !                                                                      !
                !!                                                                      !!
                !!  A mesh for the theta and phi dimensions is established in Poseidon  !!
                !!  initialization and is here used.  Here quadrature points are mapped !!
                !!  into the space of theta element, te.                                !!
                !!                                                                      !!
                 !                                                                      !
                T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), T_xlocs)
                deltat = tlocs(te+1)-tlocs(te)



                 !                                                                      !
                !!                                                                      !!
                !!  Here the second part, the Associated Legendre Polynomial, of the    !!
                !!  spherical harmonic function is precomputed at each of the           !!
                !!  theta quadrature points, and will be called into play later.        !!
                !!                                                                      !!
                 !                                                                      !
                SphereHarm_ThetaPart = SphereHarm_NormFactor*Legendre_Poly(l, -m, Source_Degrees(2), T_locs)



                 !                              !
                !!                              !!
                !!  Loop over Phi Elements      !!
                !!                              !!
                 !                              !
                DO pe = 0,NUM_P_ELEMENTS - 1



                     !                                                                          !
                    !!                                                                          !!
                    !!  Here the quadrature points are mapped from (-1,1) space to phi space.   !!
                    !!                                                                          !!
                     !                                                                          !
                    P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), P_xlocs)
                    deltap = plocs(pe+1) - plocs(pe)




                     !                                      !
                    !!                                      !!
                    !!  Loop over Radial Quadrature Points  !!
                    !!                                      !!
                     !                                      !
                    DO rd = 1,Source_Degrees(1)


                          !                                                                             !
                         !!                                                                             !!
                        !!!     Precompute Lagrange Polynomial values for current quadrature point,     !!!
                        !!!     as well as precompute the radial contribution to the triple integral.   !!!
                         !!                                                                             !!
                          !                                                                             !
                        LagP = Lagrange_Poly(R_xlocs(rd), DEGREE, Poly_xlocs)
                        R_SQR = R_locs(rd)*R_locs(rd)*deltar/2.0_idp * R_weights(rd)



                         !                                      !
                        !!                                      !!
                        !!  Loop over Phi Quadrature Points     !!
                        !!                                      !!
                         !                                      !
                        DO pd = 1, Source_Degrees(3)



                              !                                                                             !
                             !!                                                                             !!
                            !!!     Precompute phi contribution to the spherical harmonic, e^(-im*phi),     !!!
                            !!!     and phi's contribution to the triple integral.                          !!!
                             !!                                                                             !!
                              !                                                                             !
                            P_PRE = CDEXP(complex(0.0_idp,-m * P_locs(pd))) * deltap/2.0_idp * P_weights(pd)




                             !                                      !
                            !!                                      !!
                            !!  Loop over Theta Quadrature Points   !!
                            !!                                      !!
                             !                                      !
                            DO td = 1, Source_Degrees(2)


                                  !                                                                             !
                                 !!                                                                             !!
                                !!!     Precompute theta contribution to the triple integral, calling the       !!!
                                !!!     already computed theta part of the spherical harmonic.                  !!!
                                 !!                                                                             !!
                                  !                                                                             !
                                T_PRE = sin(T_locs(td))* SphereHarm_ThetaPart(td) * deltat/2.0_idp * T_weights(td)





                                 !                                                  !
                                !!                                                  !!
                                !!  Loop over the Lagrange Points/Element Nodes     !!
                                !!                                                  !!
                                 !                                                  !
                                DO p = 0,DEGREE



                                      !                                                                             !
                                     !!                                                                             !!
                                    !!!     Put together all precomputed pieces of the triple integral, along with  !!!
                                    !!!     the value of the source function at the current quadarature location.   !!!
                                     !!                                                                             !!
                                      !                                                                             !
                                    Src_Array(re*DEGREE + p, m, l) = Src_Array(re*DEGREE + p, m, l)                     &
                                                                        - Source_Term_Coefficients(re,te,pe,rd,td,pd)   &
                                                                        * 4.0_idp * pi * R_SQR*LagP(p)                  &
                                                                        * T_PRE                                         &
                                                                        * P_PRE


                                END DO

                            END DO

                        END DO

                    END DO

                END DO

            END DO

        END DO

    END DO

END DO



END SUBROUTINE Triple_Integral




END MODULE Integrator_3D_Module
