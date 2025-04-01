!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
!This is the grid setup routine below. 
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



MODULE grid_SFT

 USE variables
 USE flows
 
IMPLICIT NONE

CONTAINS

 SUBROUTINE ReadfromUser(parameterFile)
 CHARACTER*(*), INTENT(IN):: parameterFile
 

 NAMELIST /user/ dataDir, input_files, nthUnif, nphUnif,  &
                 L, eta, tau, stopDay, C, peak_lat, total_bipoles, bipolefile, &
                 savesources,ADDBIPOLES,writefluximbalance, &
                 mc1,saverestart,restartfreq,restartDir, &
                 restart, restartDay, stop_res

 OPEN(111, FILE=TRIM(parameterFile),STATUS='old')
 READ(111,nml=user) 
 CLOSE(111)
 
 CALL System("mkdir -p "//TRIM(dataDir))
 bmrdir = TRIM(dataDir)//TRIM('/BMRs')
 IF (savesources) THEN
 CALL System("mkdir -p "//TRIM(bmrdir))
 END IF
 IF (saverestart) THEN
 CALL System("mkdir -p "//TRIM(restartDir))
 END IF
 
 END SUBROUTINE ReadfromUser

SUBROUTINE setup_grid
 ALLOCATE(sg(0:nthUnif-2))
 ALLOCATE(sc(0:nthUnif-1))
 ALLOCATE(phc(0:nphUnif-1))

 ALLOCATE(sg1(1:nthUnif-1))
 ALLOCATE(sc1(1:nthUnif))
 ALLOCATE(phc1(1:nphUnif))

 ALLOCATE(MC_vel(0:nthUnif-2))

ds = 2.0_dp/nthUnif
dphi = 2*pi/nphUnif

year = 365.25_dp
day = 86400.0_dp

CALL linspace1(-1 + ds, 1 - ds, sg1)
CALL linspace1(-1 + 0.5*ds, 1 - 0.5*ds, sc1) ! excluding ghost cells
CALL linspace1(0.5*dphi, 2*pi - 0.5*dphi, phc1)

sc(0:nthUnif-1) = sc1
sg(0:nthUnif-2) = sg1
phc(0:nphUnif-1) = phc1

END SUBROUTINE setup_grid


SUBROUTINE linspace1(a1, a2, array1)
REAL(dp), INTENT(in) :: a1, a2
REAL(dp), INTENT(out) :: array1(:)
REAL(dp) :: range2
INTEGER :: n, i
n = size(array1,1)
range2 = a2 - a1
array1(1) = a1


do i=1, n
    array1(i) = a1 + range2 * (i - 1) / (n - 1)
end do
end subroutine linspace1


END MODULE grid_SFT