   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Vector_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the routines and functions associated with the source vectors for  !##!
!##!        the linear systems.  This includes the functions which perform the      !##!
!##!        integrals required to calcuate the source vector.                       !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Source_Vector                                              !##!
!##!    +102+   Deallocate_Source_Vector                                            !##!
!##!    +103+   Generate_Source_Vector                                              !##!
!##!                                                                                !##!
!##!    +201+   Triple_Integral                                                     !##!
!##!    +202+   MacLaurin_Triple_Integral                                           !##!
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
            ONLY :  idp, pi, twopi


USE Global_Variables_And_Parameters, &
            ONLY :  NUM_R_ELEMENTS,         &
                    NUM_P_ELEMENTS,         &
                    NUM_T_ELEMENTS,         &
                    rlocs, plocs, tlocs,    &
                    Source_Function_Flag,   &
                    NUM_R_NODES,            &
                    L_LIMIT,                &
                    LM_Length,              &
                    DEGREE,                 &
                    Source_Vector,          &
                    Source_Function,        &
                    Discontinuity_Flag,     &
                    Source_Degrees,         &
                    Source_Term_Coefficients

USE Additional_Functions_Module, &
            ONLY :  Legendre_Poly, Lagrange_Poly,   &
                    Norm_Factor, POWER, Map_From_X_Space, Map_To_X_Space

USE Functions_Quadrature, &
            ONLY :  Initialize_LG_Quadrature, Initialize_LGL_Quadrature


USE Integrator_3D_Module,   &
            ONLY :   Triple_Integral


USE Integrator_3D_MacLaurin_Module,   &
            ONLY :   MacLaurin_Triple_Integral

IMPLICIT NONE






!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS

!+101+###########################################################################!
!                                                                                !
!               Allocate_Source_Vector - Allocate the Source vector              !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Source_Vector

ALLOCATE(Source_Vector(0:NUM_R_NODES-1, 1:LM_Length))

Source_Vector=0.0_idp


END SUBROUTINE Allocate_Source_Vector




!+102+###########################################################################!
!                                                                                !
!             Deallocate_Source_Vector - Allocate the Source vector              !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Source_Vector

DEALLOCATE(Source_Vector)

END SUBROUTINE Deallocate_Source_Vector









!+103+##################################################################!
!                                                                       !
!          Generate_Source_Vector - Calculate the Source vector         !
!                                                                       !

!#######################################################################!
SUBROUTINE Generate_Source_Vector()



CALL Triple_Integral(Source_Vector)


END SUBROUTINE Generate_Source_Vector







END MODULE Source_Vector_Module
