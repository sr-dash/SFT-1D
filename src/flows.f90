!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
!This is the meridional circulation routine below. 
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

MODULE flows

 USE variables

CONTAINS

!****************************************************************
 SUBROUTINE MC_flow(C, sg, MC_vel)
 ! This populates the meridional circulation variable 
 ! with prescribed profile for the flow speed.
 REAL(dp), INTENT(in) :: C, sg(0:nthUnif-2)
 REAL(dp) :: p, du
 REAL(dp), INTENT(out) :: MC_vel(0:nthUnif-2)
  
 P = 2.33_dp
 du = C*(1E-3)*((1.0_dp+P)**(0.5_dp*(P+1)))/(P**(P*0.5_dp))
 MC_vel = du*sg*(SQRT(1.0_dp - sg**2))**(P)
 
 END SUBROUTINE MC_flow

 !****************************************************************
 SUBROUTINE MC_flow_v2(C,peak_lat, sg, MC_vel)
 ! This populates the meridional circulation variable 
 ! C is the peak flow speed from the user namelist file.
 ! peak_lat is the peak latitude of the meridional flow profile in degrees.
 ! with prescribed profile for the flow speed.
 REAL(dp), INTENT(in) :: C, sg(0:nthUnif-2)
 REAL(dp) :: P, du, peak_lat
 REAL(dp), INTENT(out) :: MC_vel(0:nthUnif-2)

 P = (1.0_dp/(DSIN(peak_lat*PI/180.0_dp)**2)) - 1.0_dp
 du = C*(1E-3)*((1.0_dp+P)**(0.5_dp*(P+1)))/(P**(P*0.5_dp))
 MC_vel = du*sg*(SQRT(1.0_dp - sg**2))**(P)
 
 END SUBROUTINE MC_flow_v2

END MODULE flows
