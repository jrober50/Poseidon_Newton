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

REAL(idp), PUBLIC   ::  Timer_Core_SourceInput
REAL(idp), PUBLIC   ::  Timer_Core_PrintResults

REAL(idp), PUBLIC   ::  Timer_Core_Initialization
REAL(idp), PUBLIC   ::  Timer_Matrix_Construction

REAL(idp), PUBLIC   ::  Timer_SrcVec_Construction
REAL(idp), PUBLIC   ::  Timer_SrcVec_SubParts
REAL(idp), PUBLIC   ::  Timer_SrcVec_Main

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

Timer_Core_SourceInput       = 0.0_idp
Timer_Core_PrintResults       = 0.0_idp

Timer_Core_Initialization       = 0.0_idp
Timer_Matrix_Construction       = 0.0_idp

Timer_SrcVec_Construction       = 0.0_idp
Timer_SrcVec_SubParts           = 0.0_idp
Timer_SrcVec_Main               = 0.0_idp

Timer_LinSlv_Total              = 0.0_idp

CALL TimerStart( Timer_Total )

END SUBROUTINE Init_Timers





!+302+##########################################################################!
!                                                                               !
!           Print_Timers                                                        !
!                                                                               !
!###############################################################################!
SUBROUTINE Finalize_Timers()

REAL(idp)               :: Total_Time

CALL TimerStop( Timer_Total )


Total_Time = Timer_Total-Timer_Initialize_Test_Problem

WRITE(*,*)
WRITE(*,'(10X,A)')'    Timers    '
WRITE(*,'(10X,A)')'--------------'
WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Timer_Total                      :', Total_Time, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Source Input                     :', Timer_Core_SourceInput, ' s'
WRITE(*,*)

!
! Initiliazation
!   + Matrix Construction
!
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Poseidon Initialization          :', Timer_Core_Initialization, ' s     '
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Matrix Construction           :', Timer_Matrix_Construction, ' s     '
WRITE(*,*)
!
! Source Vector Construction
!   +SubParts
!   +Main
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Source Vector Construction       :', Timer_SrcVec_Construction, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Vector Subparts        :', Timer_SrcVec_SubParts, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Vector Main            :', Timer_SrcVec_Main, ' s'
WRITE(*,*)

! Linear Solver
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Linear Solve                     :', Timer_LinSlv_Total, ' s'
WRITE(*,*)

WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Print Results                    :', Timer_Core_PrintResults, ' s'
WRITE(*,*)
WRITE(*,*)
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Missing Time                     :'  &
        ,Total_Time                     &
        - Timer_Core_Initialization     &
        - Timer_Core_SourceInput        &
        - Timer_SrcVec_Construction     &
        - Timer_LinSlv_Total            &
        - Timer_Core_PrintResults

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
