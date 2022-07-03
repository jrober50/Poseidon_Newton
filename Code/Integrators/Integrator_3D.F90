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
                    LM_Length,                  &
                    DEGREE,                     &
                    Source_Vector,              &
                    Source_Function,            &
                    Discontinuity_Flag,         &
                    Source_Degrees,             &
                    Source_Term_Coefficients,   &
                    Source_Terms

USE Additional_Functions_Module, &
            ONLY :  Legendre_Poly, Lagrange_Poly,   &
                    Norm_Factor, POWER, Map_From_X_Space, Map_To_X_Space
                    
USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature,       &
                    Initialize_LGL_Quadrature,      &
                    Initialize_Trapezoid_Quadrature

USE Timers_Module, &
            ONLY :  TimerStart,                     &
                    TimerStop,                      &
                    Timer_SourceVector_Subparts,          &
                    Timer_SourceVector_Main

USE OMP_LIB

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


COMPLEX(idp), DIMENSION(0:NUM_R_NODES-1, 1:LM_Length), INTENT(INOUT)     ::  Src_Array



COMPLEX(idp)                                        ::  Tmp_Val


INTEGER                                             ::  lm,l,m, re, rd, te, td, pe, pd, p


REAL(idp)                                           ::  drot, dtot, dpot
REAL(idp)                                           ::  SphereHarm_NormFactor

REAL(idp), ALLOCATABLE, DIMENSION(:)                ::  R_xlocs, R_locs, R_weights, &
                                                        P_xlocs, P_locs, P_weights, &
                                                        T_xlocs, T_locs, T_weights
                                                        


REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_xlocs
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  Poly_weights
REAL(KIND = idp), DIMENSION(0:DEGREE)               ::  LagP

INTEGER                                             ::  There, Here


REAL(idp),           DIMENSION(:,:),   ALLOCATABLE  ::  R_Pre
REAL(idp),           DIMENSION(:,:,:), ALLOCATABLE  ::  T_Pre
COMPLEX(KIND = idp), DIMENSION(:,:,:), ALLOCATABLE  ::  P_Pre

COMPLEX(idp), DIMENSION(1:Source_Degrees(1)*Source_Degrees(2)*Source_Degrees(3))        :: Int_Term


ALLOCATE(R_xlocs(1:Source_Degrees(1) ) )
ALLOCATE(T_xlocs(1:Source_Degrees(2) ) )
ALLOCATE(P_xlocs(1:Source_Degrees(3) ) )

ALLOCATE(R_locs(1:Source_Degrees(1) ) )
ALLOCATE(T_locs(1:Source_Degrees(2) ) )
ALLOCATE(P_locs(1:Source_Degrees(3) ) )

ALLOCATE(R_weights(1:Source_Degrees(1) ) )
ALLOCATE(T_weights(1:Source_Degrees(2) ) )
ALLOCATE(P_weights(1:Source_Degrees(3) ) )


ALLOCATE( R_Pre(1:Source_Degrees(1),0:Num_R_Nodes-1) )
ALLOCATE( T_Pre(1:Source_Degrees(2),0:Num_T_Elements-1,1:LM_Length) )
ALLOCATE( P_Pre(1:Source_Degrees(3),0:Num_P_Elements-1,1:LM_Length) )




CALL Initialize_LGL_Quadrature(DEGREE, Poly_xlocs, Poly_weights)

CALL Initialize_LG_Quadrature(Source_Degrees(1), R_xlocs, R_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(2), T_xlocs, T_weights)
CALL Initialize_LG_Quadrature(Source_Degrees(3), P_xlocs, P_weights)
!CALL Initialize_Trapezoid_Quadrature(Source_Degrees(3), P_xlocs, P_weights)

P_locs = Map_From_X_Space(0.0_idp, 2*pi, P_xlocs)
T_locs = Map_From_X_Space(0.0_idp, pi, T_xlocs)


#if defined(POSEIDON_OPENMP_FLAG)
!$OMP PARALLEL
PRINT*,"Hello From Process ",OMP_GET_THREAD_NUM()
!$OMP END PARALLEL
#endif

#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP TARGET ENTER DATA &
    !$OMP MAP( to:     ) &
    !$OMP MAP( alloc: iErr )
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC ENTER DATA &
    !$ACC COPYIN(      ) &
    !$ACC CREATE(     iErr )
#endif

CALL TimerStart( Timer_SourceVector_Subparts )

DO re = 0, Num_R_Elements -1
DO  p = 0, Degree

    There = re*Degree + p
    R_locs = Map_From_X_Space(rlocs(re), rlocs(re+1), R_xlocs)
    drot = 0.5_idp *(rlocs(re+1) - rlocs(re))

    DO rd = 1, Source_Degrees(1)

        LagP = Lagrange_Poly(R_xlocs(rd), DEGREE, Poly_xlocs)

        R_Pre(rd,There) = 4.0_idp * pi        &
                        * R_locs(rd)        &
                        * R_locs(rd)        &
                        * drot              &
                        * R_weights(rd)     &
                        * LagP(p)

    END DO
END DO
END DO




DO l = 0, L_LIMIT
DO m = -l,l
DO te = 0, NUM_T_ELEMENTS - 1

    dtot = 0.5_idp *(tlocs(te+1) - tlocs(te))
    T_locs = Map_From_X_Space(tlocs(te), tlocs(te+1), T_xlocs)
    SphereHarm_NormFactor = POWER(-1.0_idp, m)*Norm_Factor(l,-m)
    lm = l*(l+1)+m+1

    T_Pre(:,te,lm) = sin(T_locs(:))                                     &
                    * SphereHarm_NormFactor                             &
                    * Legendre_Poly(l, -m, Source_Degrees(2), T_locs)   &
                    * dtot                                              &
                    * T_weights(:)

!    DO td = 1,Source_Degrees(2)
!        PRINT*,td,te,lm,T_PRE(td,te,lm)
!    END DO

END DO ! te
END DO ! lm
END DO





DO l = 0, L_LIMIT
DO m = -l,l
DO pe = 0, Num_P_Elements - 1
DO pd = 1,Source_Degrees(3)
    lm = l*(l+1)+m+1
    dpot = 0.5_idp *(plocs(pe+1) - plocs(pe))
    P_locs = Map_From_X_Space(plocs(pe), plocs(pe+1), P_xlocs)
    P_PRE(pd,pe,lm) = CDEXP(CMPLX(0.0_idp,-m * P_locs(pd),idp)) * dpot * P_weights(pd)
END DO
END DO
END DO
END DO



CALL TimerStop( Timer_SourceVector_Subparts)
CALL TimerStart( Timer_SourceVector_Main )






Src_Array = 0.0_idp
Tmp_Val = 0.0_idp
#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(3) &
    !$OMP PRIVATE(  ) &
    !$OMP REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(3) &
    !$ACC PRIVATE(  ) &
    !$ACC PRESENT( ) &
    !$ACC REDUCTION( MIN: TimeStep )
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP PARALLEL DO SIMD COLLAPSE(3) REDUCTION(+:Src_Array) &
    !$OMP PRIVATE(  re,te,pe,lm,p,              &
    !$OMP           Here, There, TMP_Val, Int_Term   )
#endif
DO lm = 1, LM_Length
DO re = 0,NUM_R_ELEMENTS - 1
DO p = 0,DEGREE

There = re*Degree + p

DO pe = 0,NUM_P_ELEMENTS - 1
DO te = 0,NUM_T_ELEMENTS - 1

DO pd = 1, Source_Degrees(3)
DO td = 1, Source_Degrees(2)
DO rd = 1, Source_Degrees(1)


    Here = (pd-1)*Source_Degrees(2)*Source_Degrees(1)   &
         + (td-1)*Source_Degrees(1)                     &
         +  rd-1

    
!    Int_Term(Here) =  R_PRE(rd,There)               &
!                    * T_PRE(td,te,lm)               &
!                    * P_PRE(pd,pe,lm)

    Tmp_Val = Tmp_Val         &
            - Source_Terms(Here,re,te,pe)   &
            * R_PRE(rd,There)               &
            * T_PRE(td,te,lm)               &
            * P_PRE(pd,pe,lm)


END DO ! pd Loop
END DO ! td Loop
END DO ! rd Loop

!Tmp_Val = Tmp_Val -SUM( Source_Terms(:,re,te,pe)*Int_Term(:) )



END DO  ! te Loop
END DO  ! pe Loop


Src_Array(There,lm) = Src_Array(There,lm) + Tmp_Val

TMP_Val = 0.0_idp

END DO ! p Loop
END DO  ! re Loop
END DO  ! lm Loop




#if defined(POSEIDON_OPENMP_OL_FLAG)
    !$OMP END PARALLEL DO SIMD
#elif defined(POSEIDON_OPENACC_FLAG)
    !$ACC END PARALLEL LOOP
#elif defined(POSEIDON_OPENMP_FLAG)
    !$OMP END PARALLEL DO SIMD
#endif

CALL TimerStop( Timer_SourceVector_Main )


END SUBROUTINE Triple_Integral




END MODULE Integrator_3D_Module
