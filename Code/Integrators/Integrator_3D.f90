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


COMPLEX(idp), DIMENSION(0:NUM_R_NODES-1, -L_LIMIT:L_LIMIT, 0:L_LIMIT), INTENT(INOUT)     ::  Src_Array



INTEGER                                             ::  l, m, re, rd, te, td, pe, pd, d


REAL(idp)                                           ::  drot, dpot, dtot,     &
                                                        SphereHarm_NormFactor

REAL(idp), ALLOCATABLE, DIMENSION(:)                ::  R_xlocs, R_locs, R_weights, &
                                                        P_xlocs, P_locs, P_weights, &
                                                        T_xlocs, T_locs, T_weights, &
                                                        SphereHarm_ThetaPart


REAL(idp), DIMENSION(0:DEGREE)                      ::  Poly_xlocs, Poly_weights
REAL(idp), DIMENSION(1:Source_Degrees(1),0:DEGREE)  ::  LagP

REAL(idp), DIMENSION(1:Source_Degrees(1))           ::  R_SQR
REAL(idp), DIMENSION(1:Source_Degrees(2))           ::  T_PRE
COMPLEX(idp), DIMENSION(1:Source_Degrees(3))        ::  P_PRE
INTEGER                                             ::  Here
INTEGER                                             ::  There

COMPLEX(idp), DIMENSION(:), ALLOCATABLE             ::  Int_Terms

ALLOCATE( Int_Terms(1:Source_Degrees(1)*Source_Degrees(2)*Source_Degrees(3) ) )


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

dpot = 2.0_idp * pi
dtot = pi

Src_Array = 0.0_idp



DO l = 0, L_LIMIT
DO m = -l, l
DO re = 0,NUM_R_ELEMENTS - 1
DO d =  0,Degree

DO te = 0,NUM_T_ELEMENTS - 1
DO pe = 0,NUM_P_ELEMENTS - 1

    There = re*DEGREE + d

    R_locs = Map_From_X_Space(rlocs(re), rlocs(re+1), R_xlocs)
    T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), T_xlocs)
    P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), P_xlocs)


    drot = 0.50_idp*(rlocs(re+1) - rlocs(re))
    dtot = 0.50_idp*(tlocs(te+1) - tlocs(te))
    dpot = 0.50_idp*(plocs(pe+1) - plocs(pe))


    SphereHarm_NormFactor = POWER(-1.0_idp, m)*Norm_Factor(l,-m)
    SphereHarm_ThetaPart = SphereHarm_NormFactor*Legendre_Poly(l, -m, Source_Degrees(2), T_locs)



    DO rd = 1, Source_Degrees(1)
        LagP(rd,:) = Lagrange_Poly(R_xlocs(rd), DEGREE, Poly_xlocs)
    END DO

    DO rd = 1, Source_Degrees(1)
    DO td = 1, Source_Degrees(2)
    DO pd = 1, Source_Degrees(3)

        Here = (pd-1)*Source_Degrees(2)*Source_Degrees(1)       &
             + (td-1)*Source_Degrees(1)                         &
             + rd

        Int_Terms(Here) = 4.0_idp * pi *  R_locs(rd)* R_locs(rd) * drot * R_weights(rd)     &
                        * sin(T_locs(td))* SphereHarm_ThetaPart(td) * dtot* T_weights(td)   &
                        * CDEXP(complex(0.0_idp,-m * P_locs(pd))) * dpot * P_weights(pd)    &
                        * LagP(rd,d)
    END DO ! pd Loop
    END DO ! td Loop
    END DO ! rd Loop


    Src_Array(There,m,l) = -SUM(Source_Terms(:,re,te,pe)*Int_Terms(:) )




END DO  ! pe Loop
END DO  ! te Loop
END DO  ! d Loop
END DO  ! re Loop
END DO  ! m Loop
END DO  ! l Loop





END SUBROUTINE Triple_Integral




END MODULE Integrator_3D_Module
