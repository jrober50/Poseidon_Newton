   !################################################################################!
 !/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\!
!######################################################################################!
!##!                                                                                !##!
!##!                                                                                !##!
MODULE Units_Module                                                                 !##!
!##!                                                                                !##!
!##!________________________________________________________________________________!##!
!##!                                                                                !##!
!##!    Contains the functions, and subroutines used to define and convert to and   !##!
!##!    from different unit systems to that which Poseidon uses.                    !##!
!##!                                                                                !##!
!##!================================================================================!##!
!##!                                                                                !##!
!##!    Contains:                                                                   !##!
!##!                                                                                !##!
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







IMPLICIT NONE







                        !*F&S*==========================================!
                        !                                               !
                        !           Functions & Subroutines             !
                        !                                               !
                        !===============================================!

SUBROUTINE Set_Units( Units_Flag )


CHARACTER(LEN = 4)                          :: Units_Flag



IF ( Units_Flag == "C" ) THEN
!   CGS (centimeter-gram-second ) system
!
!   Mass    ->  g
!   Length  ->  cm
!   time    ->  s
!
!   c = 2.99792458e10 cm / s
!   G = 6.67428e-8    cm^3 / ( g s^2 )


    Meter = 100     ! centimeters

    Gram = 1        ! gram
    Kilogram = 1000 ! grams


    !Potential_Units = "cm/s^2"

    Grav_Constant_G = 6.67408*10**(-11) * (Meter*Meter*Meter)/(Kilogram * Second * Second)

    Speed_of_Light = 2.99792458*10**(8) * (Meter / Second )




ELSE IF ( Units_Flag == "S" ) THEN
    !   SI Units / MKS
    !
    !   Mass    ->  kg
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 2.99792458e8    m / s
    !   G = 6.67408e-11     m^2 / ( kg s^2 )


    Meter = 1       ! meter

    Gram = .001     ! kilograms
    Kilogram = 1    ! kilogram


    !Potential_Units = "m/s^2"

    Grav_Constant_G = 6.67408*10**(-11) * (Meter*Meter*Meter)/(Kilogram * Second * Second)

    Speed_of_Light = 2.99792458*10**(8) * (Meter / Second )




ELSE IF ( Units_Flag == "G" ) THEN
    !   Geometrized Units -> Poseidon Units
    !
    !   Mass    ->  g
    !   Length  ->  m
    !   time    ->  s
    !
    !   c = 1
    !   G = 1

    Meter = 1       ! meter

    Gram = 1        ! gram
    Kilogram = 1000 ! grams


    Grav_Constant_G = 1
    Speed_of_Light = 1


END IF






END SUBROUTINE Set_Units




END MODULE Units_Module