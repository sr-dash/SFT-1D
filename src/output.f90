!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
!This is the routine to write text files (if needed). 
!For detailed theory of the setup refer to the doc file.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023

! Copyright (C) Soumyaranjan Dash, University of Hawaii

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.




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
