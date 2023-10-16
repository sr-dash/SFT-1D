!<---------------------------------------------------------------------->!
!				Surface Flux Transport Model 1D
!
!This is the SERIAL version of the model.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023



MODULE output

 USE variables

CONTAINS

SUBROUTINE File_name(filename,num,FO)

INTEGER,INTENT(IN),OPTIONAL :: FO
INTEGER,INTENT(IN) :: num
CHARACTER*(*),INTENT(IN) :: filename
!! filename takes the Filename[character] (Please name according to the variable name.)
!! and the num (integer) is the file unit.

IF (PRESENT(FO)) THEN
OPEN(UNIT=num,FILE=TRIM(dataDir)//'/'//filename//'.dat',STATUS='UNKNOWN')
END IF

END SUBROUTINE File_name

END MODULE output
