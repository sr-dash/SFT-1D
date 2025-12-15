MODULE write_data
  USE netcdf

  USE variables
  USE grid_SFT

  IMPLICIT NONE

CONTAINS

SUBROUTINE Check(status)
  INTEGER, INTENT(IN) :: status
  IF (status /= NF90_NOERR) THEN
    PRINT *, "NetCDF error:", nf90_strerror(status)
    STOP 2
  END IF
END SUBROUTINE Check

!---------------------------------------
SUBROUTINE DipoleNetCDF(filename, time, dipole)
  CHARACTER(*), INTENT(IN) :: filename
  REAL(dp),     INTENT(IN) :: time(:)        ! (ntime)
  REAL(dp),     INTENT(IN) :: dipole(:,:)    ! (ntime, 2)

  INTEGER :: ncid, dim_time, dim_comp
  INTEGER :: var_time, var_dipole
  INTEGER, DIMENSION(2) :: dimids

  CALL Check(nf90_create(TRIM(filename)//'.nc', NF90_CLOBBER, ncid))

  CALL Check(nf90_def_dim(ncid, "time", SIZE(time), dim_time))
  CALL Check(nf90_def_dim(ncid, "component", SIZE(dipole, 2), dim_comp))

  dimids = (/ dim_time, dim_comp /)
  CALL Check(nf90_def_var(ncid, "time",   NF90_DOUBLE, (/ dim_time /), var_time))
  CALL Check(nf90_def_var(ncid, "dipole", NF90_DOUBLE, dimids, var_dipole))

  CALL Check(nf90_enddef(ncid))
  CALL Check(nf90_put_var(ncid, var_time,   time))
  CALL Check(nf90_put_var(ncid, var_dipole, dipole))

  CALL Check(nf90_close(ncid))
END SUBROUTINE DipoleNetCDF

!---------------------------------------
SUBROUTINE BflyNetCDF(filename, latitudes, time, bfly)
  CHARACTER(*), INTENT(IN) :: filename
  REAL(dp),     INTENT(IN) :: latitudes(:)   ! (nlat)
  REAL(dp),     INTENT(IN) :: time(:)        ! (ntime)
  REAL(dp),     INTENT(IN) :: bfly(:,:)      ! (ntime, nlat)

  INTEGER :: ncid, dim_time, dim_lat
  INTEGER :: var_lat, var_time, var_bfly
  INTEGER, DIMENSION(2) :: dimids

  CALL Check(nf90_create(TRIM(filename)//'.nc', NF90_CLOBBER, ncid))

  CALL Check(nf90_def_dim(ncid, "lat",  SIZE(latitudes), dim_lat))
  CALL Check(nf90_def_dim(ncid, "time", SIZE(time),      dim_time))

  CALL Check(nf90_def_var(ncid, "lat",  NF90_DOUBLE, (/ dim_lat /), var_lat))
  CALL Check(nf90_def_var(ncid, "time", NF90_DOUBLE, (/ dim_time /), var_time))

  dimids = (/ dim_time, dim_lat /)
  CALL Check(nf90_def_var(ncid, "bfly", NF90_DOUBLE, dimids, var_bfly))

  CALL Check(nf90_enddef(ncid))

  CALL Check(nf90_put_var(ncid, var_lat,  latitudes))
  CALL Check(nf90_put_var(ncid, var_time, time))
  CALL Check(nf90_put_var(ncid, var_bfly, bfly))

  CALL Check(nf90_close(ncid))
END SUBROUTINE BflyNetCDF

!---------------------------------------
SUBROUTINE BrMapNetCDF(filename, latitudes, br_map)
  CHARACTER(*), INTENT(IN) :: filename
  REAL(dp),     INTENT(IN) :: latitudes(:)     ! (nlat)
  REAL(dp),     INTENT(IN) :: br_map(:,:)      ! (ntime, nlat)

  INTEGER :: ncid, dim_time, dim_lat
  INTEGER :: var_lat, var_br
  INTEGER, DIMENSION(2) :: dimids

  CALL Check(nf90_create(TRIM(filename)//'.nc', NF90_CLOBBER, ncid))

  CALL Check(nf90_def_dim(ncid, "lat",  SIZE(latitudes), dim_lat))
  CALL Check(nf90_def_dim(ncid, "time", SIZE(br_map, 1), dim_time))

  CALL Check(nf90_def_var(ncid, "lat", NF90_DOUBLE, (/ dim_lat /), var_lat))

  dimids = (/ dim_time, dim_lat /)
  CALL Check(nf90_def_var(ncid, "Br", NF90_DOUBLE, dimids, var_br))

  CALL Check(nf90_enddef(ncid))

  CALL Check(nf90_put_var(ncid, var_lat, latitudes))
  CALL Check(nf90_put_var(ncid, var_br,  br_map))

  CALL Check(nf90_close(ncid))
END SUBROUTINE BrMapNetCDF

!---------------------------------------
SUBROUTINE BMRNetCDF(filename, sth, phi, bmr_map)
  CHARACTER(*), INTENT(IN) :: filename
  REAL(dp),     INTENT(IN) :: sth(:)         ! sin(theta), length nth
  REAL(dp),     INTENT(IN) :: phi(:)         ! longitude, length nph
  REAL(dp),     INTENT(IN), OPTIONAL :: bmr_map(:,:) ! shape (nth, nph)

  INTEGER :: ncid, dim_th, dim_ph
  INTEGER :: var_th, var_ph, var_bmr
  INTEGER, DIMENSION(2) :: dimids
  LOGICAL :: has_bmr

  has_bmr = PRESENT(bmr_map)

  CALL Check(nf90_create(TRIM(filename)//'.nc', NF90_CLOBBER, ncid))

  CALL Check(nf90_def_dim(ncid, "sth", SIZE(sth), dim_th))
  CALL Check(nf90_def_dim(ncid, "phi", SIZE(phi), dim_ph))

  CALL Check(nf90_def_var(ncid, "sth", NF90_DOUBLE, (/ dim_th /), var_th))
  CALL Check(nf90_def_var(ncid, "phi", NF90_DOUBLE, (/ dim_ph /), var_ph))

  IF (has_bmr) THEN
    dimids = (/ dim_th, dim_ph /)
    CALL Check(nf90_def_var(ncid, "bmr", NF90_DOUBLE, dimids, var_bmr))
  END IF

  CALL Check(nf90_enddef(ncid))

  CALL Check(nf90_put_var(ncid, var_th, sth))
  CALL Check(nf90_put_var(ncid, var_ph, phi))

  IF (has_bmr) THEN
    CALL Check(nf90_put_var(ncid, var_bmr, bmr_map))
  END IF

  CALL Check(nf90_close(ncid))
END SUBROUTINE BMRNetCDF


SUBROUTINE WriteRestart(dirname, step, C1_val, eta_val, br)
  USE variables, ONLY: nthUnif, dp
  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: dirname
  INTEGER, INTENT(IN) :: step
  REAL(dp), INTENT(IN) :: C1_val, eta_val
  REAL(dp), DIMENSION(0:nthUnif-1), INTENT(IN) :: br

  CHARACTER(5) :: step_str
  CHARACTER(200) :: filename_param, filename_br
  INTEGER :: i

  ! Format timestep as a string
  WRITE(step_str, '(I5.5)') step

  ! Construct file paths
  filename_param = TRIM(dirname)//'/restart_'//TRIM(step_str)//'.txt'
  filename_br = TRIM(dirname)//'/br_'//TRIM(step_str)//'.dat'

  ! Write parameter file
  OPEN(UNIT=91, FILE=filename_param, STATUS='unknown', ACTION='write')
    WRITE(91, '(A)') 'SFT Restart File'
    WRITE(91, *) 'C1 = ', C1_val
    WRITE(91, *) 'eta = ', eta_val
    WRITE(91, *) 'br_file = ', TRIM(filename_br)
  CLOSE(91)

  ! Write radial field (br_1D) to binary/text file
  OPEN(UNIT=92, FILE=filename_br, STATUS='unknown', ACTION='write')
    DO i = 0, nthUnif - 1
      WRITE(92, *) br(i)
    END DO
  CLOSE(92)

END SUBROUTINE WriteRestart


END MODULE write_data
