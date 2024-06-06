!<------------------------------------------------------------------>!
!				Surface Flux Transport model 1D
! This is the routine that defines all the variables used in the code. 
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



MODULE variables

 IMPLICIT NONE
 CHARACTER(LEN=200) :: dataDir, bipolefile, bFile
 CHARACTER(LEN=200) :: input_files, bmrdir, restartDir
 INTEGER, PARAMETER :: dp = KIND(1.D0)
 LOGICAL :: ADDBIPOLES, writefluximbalance, savesources
 LOGICAL :: mc1, saverestart, restart
 INTEGER :: restartfreq, restartDay, stop_res, stopDay
 REAL(dp), PARAMETER :: cflFact=0.1_dp
 REAL(dp), PARAMETER :: pi = 3.14159265358979323_dp
 REAL(dp), PARAMETER :: dtor = pi/180.0_dp
!REAL(dp), PARAMETER :: rtod = 180.0_dp/pi
 REAL(dp) :: C,P,du,C1
 REAL(dp), ALLOCATABLE :: br_1D(:),time_var(:)
 REAL(dp), ALLOCATABLE :: FV_flx(:)
 REAL(dp), ALLOCATABLE :: br_2D(:,:), brb(:,:),bfly(:,:)
 REAL(dp), ALLOCATABLE :: sc(:), phc(:), sg(:)
 REAL(dp), ALLOCATABLE :: sc1(:), phc1(:), sg1(:)
 REAL(dp), ALLOCATABLE :: lat0(:),lon0(:),sep0(:),tilt0(:),B0(:),sharpnum(:)
 REAL(dp), ALLOCATABLE :: t_yr(:), phase(:)
 
 REAL(dp), ALLOCATABLE :: MC_vel(:),MC_vel_inf(:)
 REAL(dp), ALLOCATABLE :: arrTx(:), arrTy(:,:)
 REAL(dp) :: ds, dphi, L, eta, eta1, dt, tau, bmr_a
 REAL(dp) :: dt_eta, dt_mf
 REAL(dp) :: year, day
 REAL(dp) :: dm_1D
 INTEGER :: nthUnif, nphUnif, ndt, total_bipoles
 INTEGER :: i, j, k, k1, nsteps
 INTEGER :: i1, j1, istart, junk1,idummy,bip_start

END MODULE variables
