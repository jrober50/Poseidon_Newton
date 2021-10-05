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
            ONLY :  NUM_R_ELEMENTS,             &
                    NUM_P_ELEMENTS,             &
                    NUM_T_ELEMENTS,             &
                    rlocs, plocs, tlocs,        &
                    Source_Function_Flag,       &
                    NUM_R_NODES,                &
                    L_LIMIT,                    &
                    DEGREE,                     &
                    Source_Vector,              &
                    Source_Function,            &
                    Discontinuity_Flag,         &
                    Source_Degrees,             &
                    Source_Term_Coefficients,   &
                    Source_Terms

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


REAL(KIND = idp)                                    ::  drot, dtot, dpot
REAL(KIND = idp)                                    ::  SphereHarm_NormFactor

REAL(KIND = idp), ALLOCATABLE, DIMENSION(:)         ::  R_xlocs, R_locs, R_weights, &
                                                        P_xlocs, P_locs, P_weights, &
                                                        T_xlocs, T_locs, T_weights, &
                                                        SphereHarm_ThetaPart


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_xlocs, Poly_weights, LagP


INTEGER                                             ::  There, Here

REAL(KIND = idp)                                    ::  R_SQR, T_PRE
COMPLEX(KIND = idp)                                 ::  P_PRE

REAL(idp),DIMENSION(1:Source_Degrees(1))            ::  R_SQR_Array
REAL(idp),DIMENSION(1:Source_Degrees(2))            ::  T_Pre_Array
COMPLEX(idp), DIMENSION(1:Source_Degrees(3))        ::  P_Pre_Array

REAL(idp),DIMENSION(0:Degree,1:Source_Degrees(1))   ::  LagP_Array

COMPLEX(idp), DIMENSION(1:Source_Degrees(1)*Source_Degrees(2)*Source_Degrees(3))        :: Int_Term


ALLOCATE(R_xlocs(1:Source_Degrees(1)), R_locs(1:Source_Degrees(1)), R_weights(1:Source_Degrees(1)))
ALLOCATE(P_xlocs(1:Source_Degrees(2)), P_locs(1:Source_Degrees(2)), P_weights(1:Source_Degrees(2)))
ALLOCATE(T_xlocs(1:Source_Degrees(3)), T_locs(1:Source_Degrees(3)), T_weights(1:Source_Degrees(3)))

!ALLOCATE(P_xlocs(0:P_Degree), P_locs(0:P_Degree), P_weights(0:P_Degree))
!ALLOCATE(T_xlocs(0:T_Degree), T_locs(0:T_Degree), T_weights(0:T_Degree))
ALLOCATE(SphereHarm_ThetaPart(0:Source_Degrees(2)))




CALL Initialize_LGL_Quadrature(DEGREE, Poly_xlocs, Poly_weights)

CALL Initialize_LG_Quadrature(Source_Degrees(1), R_xlocs, R_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(2), P_xlocs, P_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(3), T_xlocs, T_weights)



P_locs = Map_From_X_Space(0.0_idp, 2*pi, P_xlocs)
T_locs = Map_From_X_Space(0.0_idp, pi, T_xlocs)


Src_Array = 0.0_idp

DO pe = 0,NUM_P_ELEMENTS - 1
DO te = 0,NUM_T_ELEMENTS - 1
DO re = 0,NUM_R_ELEMENTS - 1

DO l = 0, L_LIMIT
DO m = -l, l
DO p = 0,DEGREE




    There = re*Degree + p

    R_locs = Map_From_X_Space(rlocs(re),rlocs(re+1), R_xlocs)
    T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), T_xlocs)
    P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), P_xlocs)


    drot = 0.5_idp *(rlocs(re+1) - rlocs(re))
    dtot = 0.5_idp *(tlocs(te+1) - tlocs(te))
    dpot = 0.5_idp *(plocs(pe+1) - plocs(pe))


    SphereHarm_NormFactor = POWER(-1.0_idp, m)*Norm_Factor(l,-m)
    SphereHarm_ThetaPart = SphereHarm_NormFactor*Legendre_Poly(l, -m, Source_Degrees(2), T_locs)



    DO rd = 1, Source_Degrees(1)
    DO td = 1, Source_Degrees(2)
    DO pd = 1, Source_Degrees(3)

        LagP = Lagrange_Poly(R_xlocs(rd), DEGREE, Poly_xlocs)

        R_SQR = 4.0_idp * pi * R_locs(rd)*R_locs(rd)*drot * R_weights(rd)
        T_PRE = sin(T_locs(td))* SphereHarm_ThetaPart(td) * dtot * T_weights(td)
        P_PRE = CDEXP(complex(0.0_idp,-m * P_locs(pd))) * dpot * P_weights(pd)


        Here = (pd-1)*Source_Degrees(2)*Source_Degrees(1)   &
             + (td-1)*Source_Degrees(1)                     &
             +  rd


        Src_Array(There, m, l) = Src_Array(There, m, l)         &
                                - Source_Terms(Here,re,te,pe)   &
                                * R_SQR                         &
                                * T_PRE                         &
                                * P_PRE                         &
                                * LagP(p)

    END DO ! pd Loop
    END DO ! td Loop
    END DO ! rd Loop

END DO  ! pe Loop
END DO  ! te Loop
END DO ! p Loop
END DO  ! re Loop
END DO  ! m Loop
END DO  ! l Loop



END SUBROUTINE Triple_Integral




END MODULE Integrator_3D_Module
