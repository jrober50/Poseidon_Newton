   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Source_Vector_Module                                                         !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Source_Vector                                              !##!
!##!    +102+   Deallocate_Source_Vector                                            !##!
!##!    +103+   Generate_Source_Vector                                              !##!
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
            ONLY :  idp


USE Global_Variables_And_Parameters, &
            ONLY :  NUM_R_NODES,            &
                    LM_Length,              &
                    Source_Vector


USE Integrator_3D_Module,   &
            ONLY :   Triple_Integral


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
