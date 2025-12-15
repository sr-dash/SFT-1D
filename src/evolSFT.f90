!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
!This is the timestep calculation routine below. 
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



MODULE evolSFT

 USE variables
 USE grid_SFT
 USE flows

CONTAINS

 SUBROUTINE timestep(eta,MC_vel,ds,cflFact,dt,ndt)
 
 ! This calculates the time step for the evolution 
 ! with the user provided diffusivity and flow profile.
 
 REAL(dp), INTENT(in) :: eta,ds,cflFact
 REAL(dp), INTENT(in) :: MC_vel(0:nthUnif-2)
 REAL(dp) :: dt_eta,dt_mf,day,year !,cr_day
 REAL(dp), INTENT(out) :: dt
 INTEGER, INTENT(out) :: ndt
 
 year = 365.25_dp
 day = 86400.0_dp
 ! Number of days in a carrington rotation for Br output.
 cr_day = 28*86400.0_dp

 dt_eta = (ds**2.0_dp)/eta
 dt_mf = MINVAL(ds/ABS(MC_vel))
 dt = cflFact*DMIN1(dt_eta, dt_mf)
 
 ! Ensure dt is not larger than a day
 IF (dt .LT. day) THEN
 ndt = int(day/dt)
 dt = day/ndt
 ELSE
 ndt = 1
 dt = day
 END IF

 ! Ensure dt is not larger than a carrington rotation
 !IF (dt .GT. cr_day) THEN
 !ndt = int(cr_day/dt)
 !dt = cr_day/ndt
 !ELSE
 !ndt = 1
 !dt = cr_day
 !END IF

 END SUBROUTINE timestep

 
END MODULE evolSFT 
 
