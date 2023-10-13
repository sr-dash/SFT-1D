!<---------------------------------------------------------------------->!
!				Surface Flux Transport model 1D
!This is the netCDF output routine below. For detailed theory of the setup 
!refer to the doc file.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023


MODULE write_data

 USE variables

CONTAINS
!****************************************************************
SUBROUTINE OutputToNetCDF(filename,B)

 USE netcdf
!---
 INTEGER,INTENT(IN),OPTIONAL :: B
 CHARACTER*(*),INTENT(IN) :: filename
 INTEGER(KIND=4) :: ncid   
 INTEGER(KIND=4) :: th_dimid
 INTEGER(KIND=4) :: th_varid
 INTEGER(KIND=4) :: br_varid
 !INTEGER(KIND=4) :: dm_varid
 INTEGER, PARAMETER :: dataType=NF90_REAL  ! set to NF90_REAL or NF90_DOUBLE

! ------------------------------------
! Open file:
CALL Check(nf90_create(TRIM(dataDir)//'/'//filename//'.nc', &
    NF90_CLOBBER, ncid))
! Make definitions:
CALL Check(nf90_def_dim(ncid,"th",nthUnif,th_dimid))
CALL Check(nf90_def_var(ncid,"th",dataType,th_dimid,th_varid))

IF (PRESENT(B)) THEN
CALL Check(nf90_def_var(ncid,"br",dataType,th_dimid,br_varid))
END IF

CALL Check(nf90_enddef(ncid))

! Output coordinate arrays to file:
CALL Check(nf90_put_var(ncid,th_varid,sc(0:nthUnif-1)))

IF (PRESENT(B)) THEN
! (1) Magnetic field
! ------------------
! Output variable to file:

arrTx = br_1D(0:nthUnif-1)
CALL Check(nf90_put_var(ncid,br_varid,arrTx)) 
END IF


! Close and write the file
CALL Check(nf90_close(ncid))

END SUBROUTINE OutputToNetCDF

!****************************************************************
SUBROUTINE BflyNetCDF(filename,bfly_file)

 USE netcdf
!---
 INTEGER,INTENT(IN),OPTIONAL :: bfly_file
 CHARACTER*(*),INTENT(IN) :: filename
 INTEGER(KIND=4) :: ncid   
 INTEGER(KIND=4) :: th_dimid
 INTEGER(KIND=4) :: th_varid
 INTEGER(KIND=4) :: bfly_varid
 INTEGER(KIND=4) :: time_dimid
 INTEGER(KIND=4) :: time_varid
 !INTEGER(KIND=4) :: dm_varid
 INTEGER, PARAMETER :: dataType=NF90_REAL  ! set to NF90_REAL or NF90_DOUBLE

! ------------------------------------
! Open file:
CALL Check(nf90_create(TRIM(dataDir)//'/'//filename//'.nc', &
    NF90_CLOBBER, ncid))
! Make definitions:
CALL Check(nf90_def_dim(ncid,"sth",nthUnif,th_dimid))
CALL Check(nf90_def_var(ncid,"sth",dataType,th_dimid,th_varid))
CALL Check(nf90_def_dim(ncid,"time",nsteps+1,time_dimid))
CALL Check(nf90_def_var(ncid,"time",dataType,time_dimid,time_varid))

IF (PRESENT(bfly_file)) THEN
CALL Check(nf90_def_var(ncid,"bfly",dataType, &
     (/ time_dimid, th_dimid /),bfly_varid))
END IF

CALL Check(nf90_enddef(ncid))

! Output coordinate arrays to file:
CALL Check(nf90_put_var(ncid,th_varid,sc(0:nthUnif-1)))
CALL Check(nf90_put_var(ncid,time_varid,time_var(0:nsteps)))

IF (PRESENT(bfly_file)) THEN
! (1) Butterfly diagram
! ------------------
! Output variable to file:

arrTy = bfly(0:nsteps,0:nthUnif-1)
CALL Check(nf90_put_var(ncid,bfly_varid,arrTy)) 
END IF


! Close and write the file
CALL Check(nf90_close(ncid))

END SUBROUTINE BflyNetCDF


!****************************************************************
SUBROUTINE BMRNetCDF(filename,bmr_file)

 USE netcdf
!---
 INTEGER,INTENT(IN),OPTIONAL :: bmr_file
 CHARACTER*(*),INTENT(IN) :: filename
 INTEGER(KIND=4) :: ncid   
 INTEGER(KIND=4) :: th_dimid
 INTEGER(KIND=4) :: th_varid
 INTEGER(KIND=4) :: bmr_varid
 INTEGER(KIND=4) :: ph_dimid
 INTEGER(KIND=4) :: ph_varid
 !INTEGER(KIND=4) :: dm_varid
 INTEGER, PARAMETER :: dataType=NF90_REAL  ! set to NF90_REAL or NF90_DOUBLE

! ------------------------------------
! Open file:
CALL Check(nf90_create(TRIM(bmrdir)//'/'//filename//'.nc', &
    NF90_CLOBBER, ncid))
! Make definitions:
CALL Check(nf90_def_dim(ncid,"sth",nthUnif,th_dimid))
CALL Check(nf90_def_var(ncid,"sth",dataType,th_dimid,th_varid))
CALL Check(nf90_def_dim(ncid,"phi",nphUnif,ph_dimid))
CALL Check(nf90_def_var(ncid,"phi",dataType,ph_dimid,ph_varid))

IF (PRESENT(bmr_file)) THEN
CALL Check(nf90_def_var(ncid,"bmr",dataType, &
     (/ th_dimid, ph_dimid /),bmr_varid))
END IF

CALL Check(nf90_enddef(ncid))

! Output coordinate arrays to file:
CALL Check(nf90_put_var(ncid,th_varid,sc(0:nthUnif-1)))
CALL Check(nf90_put_var(ncid,ph_varid,phc(0:nphUnif-1)))

IF (PRESENT(bmr_file)) THEN
! (1) Butterfly diagram
! ------------------
! Output variable to file:

arrTy = brb(0:nthUnif-1,0:nphUnif-1)
CALL Check(nf90_put_var(ncid,bmr_varid,arrTy)) 
END IF


! Close and write the file
CALL Check(nf90_close(ncid))

END SUBROUTINE BMRNetCDF

!****************************************************************
SUBROUTINE Check(istatus)
!----------------------------------------------------------------
! Check (ever so slightly modified from www.unidata.ucar.edu).
! For netcdf.
!----------------------------------------------------------------
 USE netcdf
 INTEGER, INTENT (IN) :: istatus
    
  IF (istatus /= nf90_noerr) THEN
     write(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
  END IF
  
END SUBROUTINE Check

END MODULE write_data
