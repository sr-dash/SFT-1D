!<---------------------------------------------------------------------->!
!				Surface Flux Transport model 1D
!This is the netCDF output routine below. For detailed theory of the setup 
!refer to the doc file.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023


MODULE variables

 IMPLICIT NONE
 CHARACTER(LEN=200):: dataDir,bipolefile,input_files,bmrdir
 INTEGER, PARAMETER :: dp = KIND(1.D0)
 LOGICAL :: ADDBIPOLES,writefluximbalance,savesources
 REAL(dp), PARAMETER :: cflFact=0.1_dp
 REAL(dp), PARAMETER :: pi = 3.14159265358979323_dp
 REAL(dp), PARAMETER :: dtor = pi/180.0_dp
!REAL(dp), PARAMETER :: rtod = 180.0_dp/pi
 REAL(dp) :: C,P,du
 REAL(dp), ALLOCATABLE :: br_1D(:),time_var(:)
 REAL(dp), ALLOCATABLE :: FV_flx(:)
 REAL(dp), ALLOCATABLE :: br_2D(:,:), brb(:,:),bfly(:,:)
 REAL(dp), ALLOCATABLE :: sc(:), phc(:), sg(:)
 REAL(dp), ALLOCATABLE :: sc1(:), phc1(:), sg1(:)
 REAL(dp), ALLOCATABLE :: lat0(:),lon0(:),sep0(:),tilt0(:),B0(:),sharpnum(:)
 REAL(dp), ALLOCATABLE :: t_yr(:), phase(:)
 
 REAL(dp), ALLOCATABLE :: MC_vel(:)
 REAL(dp), ALLOCATABLE :: arrTx(:), arrTy(:,:)
 REAL(dp) :: ds, dphi, L, eta, dt, tau, bmr_a
 REAL(dp) :: dt_eta, dt_mf
 REAL(dp) :: year, day
 REAL(dp) :: dm_1D
 INTEGER :: nthUnif, nphUnif, ndt, total_bipoles
 INTEGER :: i, j, k, k1, nsteps
 INTEGER :: i1, j1

END MODULE variables