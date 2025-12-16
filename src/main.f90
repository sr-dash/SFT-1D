
!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
!This is the main SFT function below. 
!For detailed theory of the setup refer to the doc file.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023

! Copyright Soumyaranjan Dash, University of Hawaii
! Licensed under the GNU General Public License v3.0 or later.
!
! This program comes with ABSOLUTELY NO WARRANTY.

PROGRAM SFT_1D

  ! Use custom modules for parameters, grid setup, I/O, flow fields, etc.
  USE variables
  USE grid_SFT
  USE write_data
  USE output
  USE flows
  USE init_condition
  USE evolSFT

  IMPLICIT NONE

  CHARACTER(64) :: parameterFile, restartedoutput, restartedoutput_stop
  CHARACTER(5)  :: snap1, snap2, snap3, snap4, snap5, snap6

  CALL get_command_argument(1, parameterFile)
  CALL ReadfromUser(parameterFile)

  CALL setup_grid
  ALLOCATE(br_2D(0:nthUnif-1, 0:nphUnif-1))
  ALLOCATE(br_1D(0:nthUnif-1))
  ALLOCATE(brb(0:nthUnif-1, 0:nphUnif-1))

  CALL read_bipolefile
  OPEN(13, FILE=TRIM(input_files)//'/CR2097.dat', STATUS='old', ACTION='read')
  DO i = 0, nthUnif - 1
    READ(13, *) br_2D(i, :)
  END DO
  CLOSE(13)
  CALL balance_flux(br_2D)

  IF (restart) THEN
    WRITE(snap5, '(I5.5)') restartDay
    WRITE(restartDir,*) TRIM('restart_'//snap5)
    OPEN(21, FILE=TRIM(restartDir)//'/res_'//TRIM(snap5)//'.txt', STATUS='unknown', ACTION='read')
      READ(21, *) junk1
      READ(21, *) C1
      READ(21, *) eta1
      READ(21, *) bFile
    CLOSE(21)

    OPEN(22, FILE=TRIM(bFile), STATUS='unknown', ACTION='read')
      READ(22, *) br_1D
    CLOSE(22)

    WRITE(restartedoutput,*) TRIM('restart_'//snap5//'/')
    CALL SYSTEM('mkdir -p '//TRIM(restartedoutput))
  ELSE
    br_1D = SUM(br_2D, 2) / SIZE(br_2D, 2)
    dm_1D = 1.5_dp * SUM(br_1D * sc * ds)
  END IF

  WRITE(snap2, '(I5.3)') INT(eta)
  snap2 = ADJUSTL(snap2)
  IF (mc1) THEN
    IF (restart) THEN
      CALL MC_flow(C1, sg, MC_vel)
    ELSE
      CALL MC_flow_v2(C, peak_lat, sg, MC_vel)
    END IF
    WRITE(snap3, '(I5.3)') NINT(C * 10)
    snap3 = ADJUSTL(snap3)
  END IF

  tau = 0.0_dp * day * year
  eta = MERGE(eta1, eta, restart) / (L**2)
  MC_vel = MC_vel / L

  CALL SYSTEM('mkdir -p restart_00000')
  CALL WriteRestart('restart_00000', 0, C, eta*(L**2), br_1D)

  CALL File_name('MC_vel'//TRIM(snap3), 16, FO=1)
  DO i = 0, nthUnif - 2
    WRITE(16, *) sg(i), MC_vel(i)
  END DO
  CLOSE(16)

  MC_vel = MC_vel * (1.0_dp - sg**2)

  CALL timestep(eta, MC_vel, ds, cflFact, dt, ndt)
  ALLOCATE(FV_flx(0:nthUnif))
  OPEN(12, FILE=TRIM(dataDir)//'/DM_'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS='unknown', ACTION='write')

  IF (writefluximbalance) THEN
    IF (ADDBIPOLES) THEN
      OPEN(19, FILE=TRIM(dataDir)//'/imbalance_bipoles'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS='unknown', ACTION='write')
    ELSE
      OPEN(19, FILE=TRIM(dataDir)//'/imbalance_nobipoles'//TRIM(snap2)//'_'//TRIM(snap3)//'.dat', STATUS='unknown', ACTION='write')
    END IF
  END IF

  nsteps = MERGE(stop_res, stopDay, restart)

  istart = MERGE(restartDay, 1, restart)
  IF (restart .AND. restartDay .EQ. 0) istart = 1
  IF (restart .OR. restartDay .NE. 0) restartDay = 0

  noutputs = (nsteps - istart + 1) / output_freq + 1
  ALLOCATE(time_var(0:noutputs))
  ALLOCATE(bfly(0:noutputs, 0:nthUnif-1))
  time_var(0) = 2010.457221081451_dp + REAL(istart - restartDay, dp)/365.25_dp
  
  bfly(0,:) = br_1D

  iout = 1
  IF (restart) THEN
    DO idummy = 1, total_bipoles
      IF (istart .GT. phase(idummy)) THEN
        bip_start = idummy + 1
      ELSE
        bip_start = 1
      END IF
    END DO
  ELSE
    bip_start = 1
  END IF

  DO i = istart, nsteps
    IF (ADDBIPOLES) THEN
      brb = 0.0_dp
      IF (i .LE. INT(phase(total_bipoles))) THEN
        DO k1 = bip_start, total_bipoles
          IF (i .EQ. INT(phase(k1))) THEN
            CALL make_bipole(sc, phc, dtor*lat0(k1), dtor*lon0(k1), dtor*sep0(k1)*0.56_dp, dtor*tilt0(k1), 1.0_dp, brb)
            brb = brb * ((B0(k1)/(ds*1E10*dphi*L**2)) / SUM(DABS(brb)))

            WRITE(snap1, '(I5.5)') INT(sharpnum(k1))

            IF (savesources) THEN
              ! Construct full filename and pass arrays to the updated subroutine
              CALL BMRNetCDF('/bmr_'//snap2//'_'//snap3//'_'//TRIM(snap1), sc, phc, brb)
            END IF

            br_1D = br_1D + SUM(brb, 2) / SIZE(brb, 2)
          END IF
        END DO
      END IF
    END IF


    DO j = 1, ndt
      FV_flx = 0.0_dp
      FV_flx(1:nthUnif-1) = eta * (1.0_dp - sg**2) * (br_1D(1:) - br_1D(:nthUnif-2)) / ds
      WHERE (MC_vel(0:nthUnif-2) > 0.0_dp)
        FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp * (1.0_dp + SIGN(1._dp, MC_vel)) * MC_vel * br_1D(:nthUnif-2)
      ELSEWHERE
        FV_flx(1:nthUnif-1) = FV_flx(1:nthUnif-1) - 0.5_dp * (1.0_dp - SIGN(1._dp, MC_vel)) * MC_vel * br_1D(1:)
      END WHERE
      br_1D = br_1D + (dt/ds)*(FV_flx(1:) - FV_flx(:nthUnif-1))
    END DO

    dm_1D = 1.5_dp * SUM(br_1D * sc * ds)
    IF (MOD(i - istart, output_freq) == 0) THEN
      bfly(iout,:) = br_1D
      time_var(iout) = time_var(0) + REAL(i,dp)/365.25_dp
      WRITE(12, *) time_var(iout), dm_1D
      iout = iout + 1
    END IF

    IF (writefluximbalance) WRITE(19, *) i, SUM(br_1D)
    IF (restart .AND. i == nsteps) THEN
      WRITE(snap6, '(I5.5)') i
      restartedoutput_stop = 'restart_'//snap6
      CALL SYSTEM('mkdir -p '//TRIM(restartedoutput_stop))
      CALL WriteRestart(restartedoutput_stop, i, C, eta*(L**2), br_1D)
    END IF
  END DO

  CLOSE(12)
  IF (writefluximbalance) CLOSE(19)

  ! Final Butterfly diagram output (cleaned)
  IF (restart) THEN
    CALL BflyNetCDF(restartedoutput//'bfly_restart_'//snap5//'_'//snap2//'_'//snap3, &
                    sc(0:nthUnif-1), time_var(0:iout-1), bfly(0:iout-1, 0:nthUnif-1))
  ELSE
    CALL BflyNetCDF(TRIM(dataDir)//'/bfly_'//TRIM(snap2)//'_'//TRIM(snap3), &
                    sc(0:nthUnif-1), time_var(0:iout-1), bfly(0:iout-1, 0:nthUnif-1))
  END IF

END PROGRAM SFT_1D
