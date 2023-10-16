!<---------------------------------------------------------------------->!
!				Surface Flux Transport Model 1D
!
!This is the SERIAL version of the model.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023

MODULE init_condition
 USE variables
 USE grid_SFT

CONTAINS


 SUBROUTINE read_bipolefile

  ALLOCATE(lat0(1:total_bipoles))
  ALLOCATE(lon0(1:total_bipoles))
  ALLOCATE(sep0(1:total_bipoles))
  ALLOCATE(tilt0(1:total_bipoles))
  ALLOCATE(B0(1:total_bipoles))
  ALLOCATE(t_yr(1:total_bipoles))
  ALLOCATE(phase(1:total_bipoles))
  ALLOCATE(sharpnum(1:total_bipoles))

  OPEN(UNIT=1,FILE=TRIM(bipolefile))
  DO i = 1,total_bipoles
    READ(1,*) lat0(i),lon0(i),sep0(i),tilt0(i),B0(i),t_yr(i),phase(i),sharpnum(i)
  END DO
  CLOSE(1)
  
 END SUBROUTINE read_bipolefile

 SUBROUTINE make_bipole(sc,phc,lat0,lon0,sep0,tilt0,B0,brb)
 INTEGER, PARAMETER :: dp = KIND(1.D0)
 REAL(dp), PARAMETER :: xithresh = 9.0
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: x,y,z
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: xb, yb, zb
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: thb, phib, xi
 REAL(dp), DIMENSION(0:nthUnif-1) :: sc
 REAL(dp), DIMENSION(0:nphUnif-1) :: phc
 REAL(dp), ALLOCATABLE :: sc2(:,:), phc2(:,:), sinc2(:,:)
 REAL(dp), intent(in) :: lat0,lon0,sep0,tilt0,B0
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: brb

 ALLOCATE(sc2(0:nthUnif-1,0:nphUnif-1))
 ALLOCATE(phc2(0:nthUnif-1,0:nphUnif-1))
 ALLOCATE(sinc2(0:nthUnif-1,0:nphUnif-1))

 DO i1 = 0, nthUnif-1
 DO j1 = 0, nphUnif-1
 sc2(i1,j1) = sc(i1)
 END DO
 END DO
 
 DO i1 = 0, nthUnif-1
 DO j1 = 0, nphUnif-1
 phc2(i1,j1) = phc(j1)
 END DO
 END DO
 
 sinc2 = DSQRT(1.0_dp - sc2*sc2)

 !Cartesian Coordinates
 x = DCOS(phc2)*sinc2
 y = DSIN(phc2)*sinc2
 z = sc2


 xb = x*DCOS(lat0)*DCOS(lon0) + y*DCOS(lat0)*DSIN(lon0) + z*DSIN(lat0)
 yb = x*(-DCOS(tilt0)*DSIN(lon0) + DSIN(tilt0)*DSIN(lat0)*DCOS(lon0)) &
      + y*(DCOS(tilt0)*DCOS(lon0) + DSIN(tilt0)*DSIN(lat0)*DSIN(lon0)) &
      - z*DSIN(tilt0)*DCOS(lat0)
 zb = x*(-DSIN(tilt0)*DSIN(lon0) - DCOS(tilt0)*DSIN(lat0)*DCOS(lon0)) &
      + y*(DSIN(tilt0)*DCOS(lon0) - DCOS(tilt0)*DSIN(lat0)*DSIN(lon0)) &
      + z*DCOS(tilt0)*DCOS(lat0)

 DO i1 = 0, nthUnif-1
 DO j1 = 0, nphUnif-1
 
 IF (zb(i1,j1) .GT. 1.0_dp) THEN
 zb(i1,j1) = 1.0_dp
 ELSE IF (zb(i1,j1) .LT. -1.0_dp) THEN
 zb(i1,j1) = -1.0_dp
 END IF
 
 END DO
 END DO

 
 thb = DACOS(zb)
 phib = DATAN2D(yb,xb)*(pi/180.0_dp)
 xi = (phib**2 + 2.0_dp*(0.5_dp*pi-thb)**2)/(sep0*sep0)
 brb = -(B0)*(phib/sep0)*DEXP(-xi)

 
 DO i1 = 0, nthUnif-1
 DO j1 = 0, nphUnif-1
 
 IF (xi(i1,j1) .GT. xithresh) THEN
 brb(i1,j1) = 0.0
 END IF
 
 END DO
 END DO
 CALL balance_flux(brb)

 END SUBROUTINE make_bipole


 SUBROUTINE balance_flux(f)
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: f
 REAL(dp) :: fp, fn, fpn
 REAL(dp), DIMENSION(0:nthUnif-1,0:nphUnif-1) :: br_dummy 
 
 fp = ABS(SUM(f, MASK= f .GE. 0.0_dp))
 fn = ABS(SUM(f, MASK= f .LE. 0.0_dp))

 fpn = 0.5_dp*(fp+fn)
 br_dummy(0:nthUnif-1,0:nphUnif-1) = f

 DO i1 = 0, nthUnif-1
 DO j1 = 0, nphUnif-1
 
 IF (br_dummy(i1,j1) .GT. 0.0_dp) THEN
 br_dummy(i1,j1) = br_dummy(i1,j1)*(fpn/fp)
 ELSE IF (br_dummy(i1,j1) .LT. 0.0_dp) THEN
 br_dummy(i1,j1) = br_dummy(i1,j1)*(fpn/fn)
 END IF
 
 END DO
 END DO
 f(0:nthUnif-1,0:nphUnif-1) = br_dummy

 END SUBROUTINE balance_flux

END MODULE init_condition
