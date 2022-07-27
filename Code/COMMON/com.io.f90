MODULE COMMON_IO

USE constants, &
                    ONLY : idp



TYPE :: FILE_DATA_TYPE

    CHARACTER(LEN = 30) ::  FILE_NAME = 'UNNAMED_FILE'
    CHARACTER(LEN = 20) ::  FILE_LOCATION = 'OUTPUT/'
    CHARACTER(LEN = 4)  ::  FILE_TYPE = '.OUT'
    INTEGER             ::  FILE_UNIT = 1
    LOGICAL             ::  FILE_STATUS = .FALSE.

END TYPE FILE_DATA_TYPE

CONTAINS
!###########################################################################!
!                                                                           !
!   OPEN_FILE - Opens file and sets status to true                          !
!                                                                           !
!###########################################################################!
SUBROUTINE OPEN_FILE(FILE)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE


CHARACTER(LEN = 54)                                         ::  filename


!!! Create filename from parts !!!
filename = TRIM(FILE%FILE_LOCATION) // TRIM(FILE%FILE_NAME) // TRIM(FILE%FILE_TYPE)




!!! Open File !!!
OPEN(unit = FILE%FILE_UNIT, file = filename)


!!! Change file status
FILE%FILE_STATUS = .TRUE.





END SUBROUTINE OPEN_FILE



















!###########################################################################!
!                                                                           !
!   CLOSE_FILE - Closes file and sets status to false                       !
!                                                                           !
!###########################################################################!
SUBROUTINE CLOSE_FILE(FILE)

TYPE(FILE_DATA_TYPE),INTENT(INOUT)                          ::  FILE



!!! Close file !!!
CLOSE(unit = FILE%FILE_UNIT)


!!! Change file status !!!
FILE%FILE_STATUS = .FALSE.




END SUBROUTINE CLOSE_FILE












END MODULE COMMON_IO