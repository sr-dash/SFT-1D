!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
! This is the routine that defines NETCDF file operations.
! For detailed theory of the setup refer to the doc file.
!
! Author: Soumyaranjan Dash
! Date: Jul 14 2023

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

MODULE write_data

 USE variables
 USE grid_SFT

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

IF (restart) THEN

    ! ------------------------------------
    ! Open file:
    CALL Check(nf90_create(filename//'.nc', &
        NF90_CLOBBER, ncid))
    ! Make definitions:
    CALL Check(nf90_def_dim(ncid,"sth",nthUnif,th_dimid))
    CALL Check(nf90_def_var(ncid,"sth",dataType,th_dimid,th_varid))
    CALL Check(nf90_def_dim(ncid,"time",nsteps-restartDay+1,time_dimid))
    CALL Check(nf90_def_var(ncid,"time",dataType,time_dimid,time_varid))

    IF (PRESENT(bfly_file)) THEN
    CALL Check(nf90_def_var(ncid,"bfly",dataType, &
         (/ time_dimid, th_dimid /),bfly_varid))
    END IF

    CALL Check(nf90_enddef(ncid))

    ! Output coordinate arrays to file:
    CALL Check(nf90_put_var(ncid,th_varid,sc(0:nthUnif-1)))
    CALL Check(nf90_put_var(ncid,time_varid,time_var(0:nsteps-restartDay)))

    IF (PRESENT(bfly_file)) THEN
    ! (1) Butterfly diagram
    ! ------------------
    ! Output variable to file:

    arrTy = bfly(0:nsteps-restartDay,0:nthUnif-1)
    CALL Check(nf90_put_var(ncid,bfly_varid,arrTy)) 
    END IF


    ! Close and write the file
    CALL Check(nf90_close(ncid))

ELSE
    ! ------------------------------------
    ! Open file:
    CALL Check(nf90_create(filename//'.nc', &
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

END IF

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
