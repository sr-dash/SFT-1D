!<---------------------------------------------------------------------->!
!				Surface Flux Transport model 1D
!This is the netCDF output routine below. For detailed theory of the setup 
!refer to the doc file.
! export PATH="$PATH:/usr/local/texlive/2023/bin/universal-darwin"
!Author: Soumyaranjan Dash
!Date: Jul 14 2023


PROGRAM SFT_1D

 USE variables
 USE grid_SFT
 USE write_data
 USE output
 USE flows
 USE init_condition
 USE evolSFT
 
 IMPLICIT NONE
 
 CHARACTER(64):: parameterFile
 CHARACTER*(5) :: snap1,snap2,snap3
 
 CALL get_command_argument(1, parameterFile)
 CALL ReadfromUser(parameterFile)
 CALL setup_grid

 WRITE (snap2, FMT='(i5.3)') int(eta)
 snap2 = ADJUSTL(snap2)
 
 IF (mc1) THEN
 CALL MC_flow(C,sg,MC_vel)
 !CALL MC_flow_assym(C,sg,MC_vel)
 !MC_vel = MC_vel*(1.0_dp - sg**2)
 ELSE IF (mc2) THEN
 OPEN(UNIT=12,FILE=TRIM(mcfile))
  DO i = 0,nthUnif-2
    READ(12,*) MC_vel(i)
  END DO
  CLOSE(12)
  MC_vel = MC_vel/1E3
  !MC_vel = MC_vel*(1.0_dp - sg**2)
 ELSE IF (inflow) THEN
 OPEN(UNIT=12,FILE=TRIM(mcfile))
  DO i = 0,nthUnif-2
    READ(12,*) MC_vel_inf(i)
  END DO
  CLOSE(12)
  MC_vel_inf = MC_vel_inf/1E3
  !MC_vel = MC_vel*(1.0_dp - sg**2)
 END IF
 
 IF (mc1) THEN
 WRITE (snap3, FMT='(i5.3)') int(MAXVAL(MC_vel)*10E3)
 snap3 = ADJUSTL(snap3)
 ! PRINT*,snap3
 ELSE IF (mc2) THEN
 snap3 = ADJUSTL('zfit')
 ELSE IF (inflow) THEN
 snap3 = ADJUSTL('inf')
 END IF
 
 tau = 0.0_dp*day*year
 eta = eta/(L**2)
 MC_vel = Mc_vel/L
 IF (inflow) THEN
 MC_vel_inf = Mc_vel_inf/L
 MC_vel_inf = MC_vel_inf*(1.0_dp - sg**2)
 END IF

 
 CALL File_name('MC_vel'//TRIM(snap3),16,FO=1)
 DO i=0,nthUnif-2
 WRITE(16,*)sg(i),MC_vel(i)
 END DO
 CLOSE(16)

 MC_vel = MC_vel*(1.0_dp - sg**2)


 CALL timestep(eta,MC_vel,ds,cflFact,dt,ndt)
 !IF (inflow) THEN
 !CALL timestep_inf(eta,MC_vel,MC_vel_inf,ds,cflFact,dt,ndt)
 !ELSE
 !CALL timestep(eta,MC_vel,ds,cflFact,dt,ndt)
 !END IF
 
ALLOCATE(br_2D(0:nthUnif-1,0:nphUnif-1))
ALLOCATE(br_1D(0:nthUnif-1))

ALLOCATE(brb(0:nthUnif-1,0:nphUnif-1))

 CALL read_bipolefile
 
 OPEN(13, FILE=TRIM(input_files)//'/CR2097.dat', STATUS="old", ACTION="read")
 DO i = 0,nthUnif-1
 READ(13,*) br_2D(i,:)
 END DO
 
 CALL balance_flux(br_2D)
 
 br_1D = SUM(br_2D, 2)/SIZE(br_2D,2)
 dm_1D = 1.5_dp*SUM(br_1D*sc*ds)

!Allocate and initialize Finite-Volume Fluxes at interior ribs (size:181)
ALLOCATE(FV_flx(0:nthUnif)) 
 OPEN(12, FILE=TRIM(dataDir)//'/DM_'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS="unknown", ACTION="write")
IF (writefluximbalance) THEN
    IF (ADDBIPOLES) THEN
 OPEN(19, FILE=TRIM(dataDir)//'/imbalance_bipoles'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS="unknown", ACTION="write")
    ELSE
 OPEN(19, FILE=TRIM(dataDir)//'/imbalance_nobipoles'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS="unknown", ACTION="write")  
    END IF    
END IF
 FV_flx = 0.0_dp
 nsteps = 14*365
 ALLOCATE(bfly(0:nsteps,0:nthUnif-1))
 bfly(0,:) = br_1D
 ALLOCATE(time_var(0:nsteps))
 time_var(0) = 2010.457221081451_dp
 k1 = 1
 
 bmr_a = 0.56_dp

 DO i=1, nsteps

 IF (ADDBIPOLES) THEN
    brb(0:nthUnif-1,0:nphUnif-1) = 0.0_dp
    IF (i .LE. int(phase(total_bipoles))) THEN
       DO k1 = 1,total_bipoles
         ! print*,i,total_bipoles,phase(2),k1
       IF (i .EQ. int(phase(k1))) THEN
          CALL make_bipole(sc,phc,dtor*lat0(k1),dtor*lon0(k1),dtor*sep0(k1)*bmr_a,dtor*tilt0(k1),1.0_dp,brb)
          !Scale the BMR with the observed magnetic field 
          ! The multiplicative factor is for converting the flux to field (L is in Km, so extra 1E10 is multiplied)
          brb = brb*((B0(k1)/(ds*1E10*dphi*L**2))/SUM(DABS(brb)))
          WRITE (snap1, FMT='(i5.5)') int(sharpnum(k1))
          snap1 = ADJUSTL(snap1)
          IF (savesources) THEN
            CALL BMRNetCDF('/bmr_'//TRIM(snap2)//'_'//TRIM(snap3)//'_'//TRIM(snap1),bmr_file=1)
          END IF
          br_1D = br_1D + SUM(brb, 2)/SIZE(brb,2)
       END IF
       
       END DO
    END IF
 
 END IF

 ! PRINT*,"Time step: ",dt

 DO j = 1, ndt
 FV_flx(0:nthUnif) = 0.0_dp
 
 ! Diffusion term
 FV_flx(1:nthUnif-1) = eta*(1.0_dp - sg**2)*(br_1D(1:) - br_1D(:nthUnif-2))/ds
 
 IF(inflow .AND. i .GE. 168 .AND. i .LE. 198) THEN
 ! Meridional flow by Up-Winding
 WHERE(MC_vel_inf(0:nthUnif-2).GT.0.0_dp)
 FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp*(1.0_dp + &
               SIGN(1._dp,MC_vel_inf(:)))*MC_vel_inf*br_1D(:nthUnif-2)
 ELSEWHERE
 FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp*(1.0_dp - &
               SIGN(1._dp,MC_vel_inf(:)))*MC_vel_inf*br_1D(1:)
 END WHERE
 ! Update with exponential decay (if choosen)
 IF (tau .GT. 0.0_dp) THEN
 br_1D = br_1D + (dt/ds)*(FV_flx(1:) - FV_flx(:nthUnif-1)) - (dt/tau)*br_1D
 ELSE
 br_1D = br_1D + (dt/ds)*(FV_flx(1:) - FV_flx(:nthUnif-1))
 END IF
 
 ELSE

 ! Meridional flow by Up-Winding
 WHERE(MC_vel(0:nthUnif-2).GT.0.0_dp)
 FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp*(1.0_dp + &
               SIGN(1._dp,MC_vel(:)))*MC_vel*br_1D(:nthUnif-2)
 ELSEWHERE
 FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp*(1.0_dp - &
               SIGN(1._dp,MC_vel(:)))*MC_vel*br_1D(1:)
 END WHERE
 ! Update with exponential decay (if choosen)
 IF (tau .GT. 0.0_dp) THEN
 br_1D = br_1D + (dt/ds)*(FV_flx(1:) - FV_flx(:nthUnif-1)) - (dt/tau)*br_1D
 ELSE
 br_1D = br_1D + (dt/ds)*(FV_flx(1:) - FV_flx(:nthUnif-1))
 END IF
 
 END IF
 
 END DO
 
 
 dm_1D = 1.5_dp*SUM((br_1D)*sc*ds)
 bfly(i,:) = br_1D
 time_var(i) = time_var(0) + i/365.25_dp
 WRITE(12,*)i/365.25_dp, dm_1D
 IF (writefluximbalance) THEN
 WRITE(19,*)i, SUM(br_1D)
 END IF

 END DO
 CLOSE(12)
 IF (writefluximbalance) THEN
 CLOSE(19)
 END IF
 CALL BflyNetCDF('/bfly_'//TRIM(snap2)//'_'//TRIM(snap3),bfly_file=1)

END PROGRAM SFT_1D
