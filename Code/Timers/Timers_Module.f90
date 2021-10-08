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
REAL(idp), PUBLIC   ::  Timer_SourceInput_PartA
REAL(idp), PUBLIC   ::  Timer_SourceInput_PartB

REAL(idp), PUBLIC   ::  Timer_Core_PrintResults

REAL(idp), PUBLIC   ::  Timer_Core_Initialization
REAL(idp), PUBLIC   ::  Timer_Initialization_Matrix

REAL(idp), PUBLIC   ::  Timer_Core_SourceVector
REAL(idp), PUBLIC   ::  Timer_SourceVector_SubParts
REAL(idp), PUBLIC   ::  Timer_SourceVector_Main

REAL(idp), PUBLIC   ::  Timer_Core_LinearSolve


CONTAINS


!+301+##########################################################################!
!                                                                               !
!           Init_Timers                                                          !
!                                                                               !
!###############################################################################!
SUBROUTINE Init_Timers

Timer_Total                     = 0.0_idp

Timer_Initialize_Test_Problem   = 0.0_idp

Timer_Core_SourceInput          = 0.0_idp
Timer_SourceInput_PartA         = 0.0_idp
Timer_SourceInput_PartB         = 0.0_idp

Timer_Core_PrintResults         = 0.0_idp

Timer_Core_Initialization       = 0.0_idp
Timer_Initialization_Matrix     = 0.0_idp

Timer_Core_SourceVector         = 0.0_idp
Timer_SourceVector_SubParts     = 0.0_idp
Timer_SourceVector_Main         = 0.0_idp

Timer_Core_LinearSolve          = 0.0_idp

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
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Input Part A           :', Timer_SourceInput_PartA, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Input Part B           :', Timer_SourceInput_PartB, ' s'
WRITE(*,*)

!
! Initiliazation
!   + Matrix Construction
!
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Poseidon Initialization          :', Timer_Core_Initialization, ' s     '
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Matrix Construction           :', Timer_Initialization_Matrix,' s     '
WRITE(*,*)
!
! Source Vector Construction
!   +SubParts
!   +Main
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Source Vector Construction       :', Timer_Core_SourceVector, ' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Vector Subparts        :', Timer_SourceVector_SubParts,' s'
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  '  -Source Vector Main            :', Timer_SourceVector_Main, ' s'
WRITE(*,*)

! Linear Solver
WRITE(*,'(7X,A,5X,ES12.6E2,A)') &
  'Linear Solve                     :', Timer_Core_LinearSolve, ' s'
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
        - Timer_Core_SourceVector       &
        - Timer_Core_LinearSolve        &
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
