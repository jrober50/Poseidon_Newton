   !##########################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!###############################################################################!
!##!                                                                         !##!
!##!                                                                         !##!
MODULE Timers_Module                                                          !##!
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
        ONLY :  idp

USE, INTRINSIC :: ISO_FORTRAN_ENV, ONLY: &
  I8 => INT64


IMPLICIT NONE

REAL(idp), PUBLIC   ::  Timer_Total

REAL(idp), PUBLIC   ::  Timer_Initialize_Test_Problem

REAL(idp), PUBLIC   ::  Timer_Matrix_Construction

REAL(idp), PUBLIC   ::  Timer_SrcVec_Construction

REAL(idp), PUBLIC   ::  Timer_LinSlv_Total


CONTAINS


!+301+##########################################################################!
!                                                                               !
!           Init_Timers                                                          !
!                                                                               !
!###############################################################################!
SUBROUTINE Init_Timers

Timer_Total                     = 0.0_idp

Timer_Initialize_Test_Problem   = 0.0_idp

Timer_Matrix_Construction       = 0.0_idp

Timer_SrcVec_Construction       = 0.0_idp

Timer_LinSlv_Total              = 0.0_idp

CALL TimerStart( Timer_Total )

END SUBROUTINE Init_Timers





!+302+##########################################################################!
!                                                                               !
!           Print_Timers                                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Finalize_Timers()

CALL TimerStop( Timer_Total )


WRITE(*,*)
WRITE(*,'(10X,A)')'    Timers    '
WRITE(*,'(10X,A)')'--------------'
WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Timer_Total                                :', Timer_Total-Timer_Initialize_Test_Problem, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Matrix Construction                        :', Timer_Matrix_Construction, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Source Vector Construction                 :', Timer_SrcVec_Construction, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Linear Solve                               :', Timer_LinSlv_Total, ' s'
WRITE(*,*)
WRITE(*,*)
WRITE(*,*)

END SUBROUTINE Finalize_Timers




!+101+##########################################################################!
!                                                                               !
!      TimersStart                                               				!
!                                                                               !
!###############################################################################!
SUBROUTINE TimerStart( Timer )

REAL(idp), INTENT(INOUT)            :: Timer

Timer = Timer - TimerWtime()

END SUBROUTINE TimerStart




!+102+##########################################################################!
!                                                                               !
!      TimersStop                                                               !
!                                                                               !
!###############################################################################!
SUBROUTINE TimerStop( Timer )

REAL(idp), INTENT(INOUT)            :: Timer

Timer = Timer + TimerWtime()

END SUBROUTINE TimerStop




!+201+##########################################################################!
!                                                                               !
!        TimersWtime                                                            !
!                                                                               !
!###############################################################################!
REAL(idp) FUNCTION TimerWtime()

INTEGER(I8)                ::  Clock_Read
INTEGER(I8)                ::  Clock_Rate
INTEGER(I8)                ::  Clock_Max

CALL SYSTEM_CLOCK( Clock_Read, Clock_Rate, Clock_Max)
TimerWtime = REAL( Clock_Read, idp )/REAL( Clock_Rate, idp )

RETURN

END FUNCTION TimerWtime










END MODULE Timers_Module
