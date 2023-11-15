!<---------------------------------------------------------------------->!
!				Surface Flux Transport Model 1D
!
!This is the SERIAL version of the model.
!
!Author: Soumyaranjan Dash
!Date: Jul 14 2023

MODULE flows

 USE variables

CONTAINS

!****************************************************************
 SUBROUTINE MC_flow(C, sg, MC_vel)
 ! This populates the meridional circulation variable 
 ! with prescribed profile for the flow speed.
 REAL(dp), INTENT(in) :: C, sg(0:nthUnif-2)
 REAL(dp) :: P, du
 REAL(dp), INTENT(out) :: MC_vel(0:nthUnif-2)
  
 P = 2.33_dp
 du = C*((1.0_dp+P)**(0.5_dp*(P+1)))/(P**(P*0.5_dp))
 MC_vel = du*sg*(SQRT(1.0_dp - sg**2))**(P)
 
 END SUBROUTINE MC_flow
 
 !****************************************************************
 SUBROUTINE MC_flow_assym(C, sg, MC_vel)
 ! This populates the meridional circulation variable 
 ! with prescribed profile for the flow speed.
 REAL(dp), INTENT(in) :: C, sg(0:nthUnif-2)
 REAL(dp) :: P, du, assym_factor
 REAL(dp), INTENT(out) :: MC_vel(0:nthUnif-2)
 
 assym_factor = -0.5_dp
 P = 2.33_dp
 du = C*((1.0_dp+P)**(0.5_dp*(P+1)))/(P**(P*0.5_dp))
 MC_vel = (1.0_dp + assym_factor*sg)*du*sg*(SQRT(1.0_dp - sg**2))**(P)
 
 END SUBROUTINE MC_flow_assym

END MODULE flows
