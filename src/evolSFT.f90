!<---------------------------------------------------------------------->!
!				Surface Flux Transport model 1D
!This is the netCDF output routine below. For detailed theory of the setup 
!refer to the doc file.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023

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
 REAL(dp) :: dt_eta,dt_mf,day,year
 REAL(dp), INTENT(out) :: dt
 INTEGER, INTENT(out) :: ndt
 
 year = 365.25_dp
 day = 86400.0_dp

 dt_eta = (ds**2.0_dp)/eta
 dt_mf = MINVAL(ds/ABS(MC_vel))
 dt = cflFact*DMIN1(dt_eta, dt_mf)
 IF (dt .LT. day) THEN
 ndt = int(day/dt)
 dt = day/ndt
 ELSE
 ndt = 1
 dt = day
 END IF

 END SUBROUTINE timestep

 
END MODULE evolSFT 
 
