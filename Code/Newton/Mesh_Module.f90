   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Mesh_Module                                                                  !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the top level subroutines needed to inialize, run, and close       !##!
!##!        Poseidon.                                                               !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
!##!    +101+   Allocate_Mesh                                                       !##!
!##!    +102+   Deallocate_Mesh                                                     !##!
!##!                                                                                !##!
!##!    +201+   Generate_Defined_Mesh                                               !##!
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
                ONLY :  idp, pi


USE Global_Variables_And_Parameters, &
                ONLY :  NUM_R_ELEMENTS, NUM_T_ELEMENTS, NUM_P_ELEMENTS,         &
                        R_OUTER, R_INNER, Source_Function_Flag, R_MAX, R_MIN,   &
                        rlocs, tlocs, plocs,                                    &
                        Discontinuity_Flag

IMPLICIT NONE




!*F&S*==========================================!
!                                               !
!           Functions & Subroutines             !
!                                               !
!===============================================!
CONTAINS


!+101+###########################################################################!
!                                                                                !
!                     ALLOCATE_MESH - Allocate the mesh arrays                   !
!                                                                                !
!################################################################################!
SUBROUTINE Allocate_Mesh

ALLOCATE(rlocs(0:NUM_R_ELEMENTS))
ALLOCATE(tlocs(0:NUM_T_ELEMENTS))
ALLOCATE(plocs(0:NUM_P_ELEMENTS))

END SUBROUTINE Allocate_Mesh


!+102+###########################################################################!
!                                                                                !
!                  DEALLOCATE_MESH - Dellocate the mesh arrays                   !
!                                                                                !
!################################################################################!
SUBROUTINE Deallocate_Mesh

DEALLOCATE(rlocs)
DEALLOCATE(tlocs)
DEALLOCATE(plocs)

END SUBROUTINE Deallocate_Mesh









!+201+##################################################################################!
!                                                                                       !
!       GENERATE_DEFINED_MESH - Generate the values for the mesh sent in using          !
!                                predefined values for the width of each element.       !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Input:  Mesh_Start - Single Real number defining inner boundary location.           !
!                                                                                       !
!           Number_of_Elements - Single Integer value defining the number of elements   !
!                                   in the mesh.                                        !
!                                                                                       !
!           Element_Width_Vector - Real valued vector of length(1:Number_of_Elements)   !
!                                       containing Real numbers defining the width of   !
!                                        each element.                                  !
!                                                                                       !
!---------------------------------------------------------------------------------------!
!                                                                                       !
!   Output: Mesh - Real Vector,length(0:Number_of_Elements) that on output contains     !
!                       values describing the element edge locations.                   !
!                                                                                       !
!#######################################################################################!
SUBROUTINE Generate_Defined_Mesh(Number_of_Elements, Mesh_Start, Element_Width_Vector, Mesh)


INTEGER, INTENT(IN)                                                 ::  Number_of_Elements
REAL(KIND = idp), INTENT(IN)                                        ::  Mesh_Start
REAL(KIND = idp), DIMENSION(1:Number_of_Elements), INTENT(IN)       ::  Element_Width_Vector

REAL(KIND = idp), DIMENSION(0:Number_of_Elements), INTENT(OUT)      ::  Mesh

INTEGER                                                             ::  i


mesh(0) = Mesh_Start
DO i = 1,Number_of_Elements


    mesh(i) = mesh(i-1) + Element_Width_Vector(i)


END DO


END SUBROUTINE Generate_Defined_Mesh














END MODULE Mesh_Module