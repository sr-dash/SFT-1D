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
 CHARACTER(64) :: restartedoutput
 CHARACTER*(5) :: snap1,snap2,snap3,snap4,snap5
 
 CALL get_command_argument(1, parameterFile)
 CALL ReadfromUser(parameterFile)
 CALL setup_grid

ALLOCATE(br_2D(0:nthUnif-1,0:nphUnif-1))
ALLOCATE(br_1D(0:nthUnif-1))

ALLOCATE(brb(0:nthUnif-1,0:nphUnif-1))

 CALL read_bipolefile
 
 OPEN(13, FILE=TRIM(input_files)//'/CR2097.dat', STATUS="old", ACTION="read")
 DO i = 0,nthUnif-1
 READ(13,*) br_2D(i,:)
 END DO
 
 CALL balance_flux(br_2D)
 

! Restart simulation functionalities.
 IF (restart) THEN
   WRITE(snap5,FMT='(i5.5)')restartDay
   OPEN(21, FILE=TRIM(restartDir)//'/res_'//TRIM(snap5)//'.txt', STATUS="unknown", ACTION="read")
   READ(21,fmt='(I5.5)')junk1
   READ(21,fmt='(F7.5)')C1
   READ(21,fmt='(F7.3)')eta1
   READ(21,fmt='(A)')bFile
   CLOSE(21)
   OPEN(22, FILE=TRIM(bFile), STATUS="unknown", ACTION="read")
   READ(22,*)br_1D
   CLOSE(22)
   PRINT*,bFile
   IF (restart) THEN
    WRITE(restartedoutput,*)TRIM('restart_'//snap5//'/')
    CALL System("mkdir -p "//TRIM(restartedoutput))
   END IF
 ELSE
   br_1D = SUM(br_2D, 2)/SIZE(br_2D,2)
   dm_1D = 1.5_dp*SUM(br_1D*sc*ds)
 END IF



 WRITE (snap2, FMT='(i5.3)') int(eta)
 snap2 = ADJUSTL(snap2)
 
 IF (mc1) THEN
   IF (restart) THEN
      CALL MC_flow(C1,sg,MC_vel)
   ELSE
      CALL MC_flow(C,sg,MC_vel)
   END IF
 !MC_vel = MC_vel*(1.0_dp - sg**2)
 END IF
 
 IF (mc1) THEN
 WRITE (snap3, FMT='(i5.3)') int(MAXVAL(MC_vel)*10E3)
 snap3 = ADJUSTL(snap3)
 END IF
 
 tau = 0.0_dp*day*year
 IF (restart) THEN
   eta = eta1/(L**2)
ELSE
   eta = eta/(L**2)
END IF
 MC_vel = Mc_vel/L
 
 CALL File_name('MC_vel'//TRIM(snap3),16,FO=1)
 DO i=0,nthUnif-2
 WRITE(16,*)sg(i),MC_vel(i)
 END DO
 CLOSE(16)

 MC_vel = MC_vel*(1.0_dp - sg**2)


 CALL timestep(eta,MC_vel,ds,cflFact,dt,ndt)

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
 IF (restart) THEN
  ALLOCATE(bfly(0:nsteps-restartDay,0:nthUnif-1))
  bfly(0,:) = br_1D
  ALLOCATE(time_var(0:nsteps-restartDay))
  time_var(0) = 2010.457221081451_dp + restartDay/365.25_dp
  PRINT*,SIZE(bfly,1),SIZE(bfly,2),nsteps
 ELSE
  ALLOCATE(bfly(0:nsteps,0:nthUnif-1))
  bfly(0,:) = br_1D
  ALLOCATE(time_var(0:nsteps))
  time_var(0) = 2010.457221081451_dp
  PRINT*,SIZE(bfly,1),SIZE(bfly,2),nsteps
 END IF
 k1 = 1
 
 bmr_a = 0.56_dp

 IF (restart) THEN
   istart = restartDay
 ELSE
   istart = 1
 END IF

IF (restart) THEN
 DO idummy=1,total_bipoles
 IF(istart .GT. phase(idummy)) THEN
 bip_start = idummy+1
 END IF
 END DO

 DO i=istart , nsteps
 IF (ADDBIPOLES) THEN
 brb(0:nthUnif-1,0:nphUnif-1) = 0.0_dp

  IF (i .LE. int(phase(total_bipoles))) THEN
 
 DO k1 = bip_start,total_bipoles
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

 DO j = 1, ndt
 FV_flx(0:nthUnif) = 0.0_dp

 ! Diffusion term
 FV_flx(1:nthUnif-1) = eta*(1.0_dp - sg**2)*(br_1D(1:) - br_1D(:nthUnif-2))/ds
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
 
 END DO
 
 dm_1D = 1.5_dp*SUM((br_1D)*sc*ds)
!PRINT*,i-istart,nsteps
 bfly(i-istart,:) = br_1D
 time_var(i-istart) = time_var(0) + (i-istart)/365.25_dp
 WRITE(12,*)i/365.25_dp, dm_1D
 IF (writefluximbalance) THEN
 WRITE(19,*)i, SUM(br_1D)
 END IF

 END DO
 CLOSE(12)

 IF (writefluximbalance) THEN
 CLOSE(19)
 END IF

! With restart = .FALSE.
ELSE

 DO i=1 , nsteps
 !PRINT*,i,nsteps
 IF (ADDBIPOLES) THEN
 brb(0:nthUnif-1,0:nphUnif-1) = 0.0_dp

    IF (i .LE. int(phase(total_bipoles))) THEN
    
       DO k1 = 1,total_bipoles
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

 DO j = 1, ndt
    FV_flx(0:nthUnif) = 0.0_dp

    ! Diffusion term
    FV_flx(1:nthUnif-1) = eta*(1.0_dp - sg**2)*(br_1D(1:) - br_1D(:nthUnif-2))/ds
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
 
 END DO
 
 IF (saverestart) THEN
   IF (MODULO(i,restartfreq) .EQ. 0) THEN
      WRITE (snap4, FMT='(i5.5)') i
      OPEN(15, FILE=TRIM(restartDir)//'/res_'//TRIM(snap4)//'.txt', STATUS="unknown", ACTION="write")
      WRITE(15,fmt='(I5.5)')i
      WRITE(15,fmt='(F7.5)')C
      WRITE(15,fmt='(F7.3)')eta*(L**2)
      WRITE(15,fmt='(A)')TRIM(TRIM(restartDir)//'/b_'//TRIM(snap4)//'.txt')
      CLOSE(15)
      OPEN(17, FILE=TRIM(restartDir)//'/b_'//TRIM(snap4)//'.txt', STATUS="unknown", ACTION="write")
      WRITE(17, *)br_1D
      CLOSE(17)
   END IF
 END IF
 
 dm_1D = 1.5_dp*SUM((br_1D)*sc*ds)
!PRINT*,i-istart,nsteps
 bfly(i,:) = br_1D
 time_var(i) = time_var(0) + (i)/365.25_dp
 WRITE(12,*)i/365.25_dp, dm_1D
 IF (writefluximbalance) THEN
 WRITE(19,*)i, SUM(br_1D)
 END IF

 END DO
 CLOSE(12)


 IF (writefluximbalance) THEN
 CLOSE(19)
 END IF


END IF


 IF (restart) THEN
   PRINT*,TRIM(restartedoutput)//TRIM('bfly_')//TRIM(snap2)//'_'//TRIM(snap3)
   CALL BflyNetCDF(TRIM(restartedoutput)//TRIM('bfly_restart_')//TRIM(snap5)//'_'//TRIM(snap2)//'_'//TRIM(snap3),bfly_file=1)
 ELSE
   PRINT*,TRIM(dataDir)//'/bfly_'//TRIM(snap2)//'_'//TRIM(snap3) 
   CALL BflyNetCDF(TRIM(dataDir)//'/bfly_'//TRIM(snap2)//'_'//TRIM(snap3),bfly_file=1)
 END IF

END PROGRAM SFT_1D
