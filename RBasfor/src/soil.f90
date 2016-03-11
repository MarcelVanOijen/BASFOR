module soil

use parameters
use management
use environment
implicit none

real :: Drain, Runoff, WAFC, WAST, WC, WCFC
real :: dCLITT, rCLITT, rCSOMF, Rsoil
real :: dCLITTrsoil, dCLITTsomf, dCSOMF, dCSOMFrsoil, dCSOMFsoms, dCSOMS
real :: Nemission, NemissionN2O, NemissionNO, Nfixation, Nleaching
real :: NLITTnmin, NLITTsomf, NlossTreeLB, NlossTreeR, Nmineralisation
real :: dNLITT, dNSOMF, dNSOMS, NSOMFnmin, NSOMFsoms, rNLITT, rNSOMF
real :: Tsoil, fTsoil, Thist(10)

Contains

  Subroutine water(WA,RAINint,Evap,Tran,LAI)
  !=============================================================================
  ! Calculate rates of runoff and drainage (mm d-1)
  ! Author - Marcel van Oijen (CEH-Edinburgh)
  ! 6-11-2005
  !=============================================================================
  integer :: day
  real    :: Evap, LAI, RAINint, Tran, WA
  WCFC   = WCST * FWCFC                                          ! % (m3 m-3)
  WC     = 0.001 * WA   / ROOTD                                  ! % (m3 m-3)
  WAFC   = 1000. * WCFC * ROOTD                                  ! % (mm)
  WAST   = 1000. * WCST * ROOTD                                  ! % (mm)
  RUNOFF = (RAIN-RAINint) * sin(atan(SLOPE/100)) * &
                                        exp(-KRUNOFF * LAI)      ! % (mm d-1)
  Drain  = max(0. , (WA-WAFC)/DELT + &
                    RAIN - RAINint - Runoff - Evap - Tran ) ! % (mm d-1)
  end Subroutine water

  Subroutine CNsoil(RWA,WFPS,WA,gCR,CLITT,CSOMF,NLITT,NSOMF,NSOMS,NMIN,CSOMS)
  real :: CLITT, CSOMF, CSOMS, fN2O, gCR, NLITT, NMIN, NSOMF, NSOMS
  real :: RWA, WA, WFPS
  
  ! C Litter
  rCLITT      = ((CLITT*Runoff) / ROOTD) * RRUNBULK * 0.001
  dCLITT      =  (CLITT*fTsoil) / TCLITT
  dCLITTsomf  = FLITTSOMF * dCLITT
  dCLITTrsoil = dCLITT - dCLITTsomf

  ! C SOM fast
  rCSOMF      = ((CSOMF*Runoff) / ROOTD) * RRUNBULK * 0.001
  dCSOMF      =  (CSOMF*fTsoil) / TCSOMF
  dCSOMFsoms  = FSOMFSOMS * dCSOMF
  dCSOMFrsoil = dCSOMF - dCSOMFSOMS
  
  ! C SOM slow
  dCSOMS      = (CSOMS*fTsoil) / TCSOMS

  ! Respiration
  Rsoil       = dCLITTrsoil + dCSOMFrsoil + dCSOMS

  ! N Litter
  rNLITT      = ((NLITT*Runoff) / ROOTD) * RRUNBULK * 0.001
  dNLITT      =  (NLITT*dCLITT) / CLITT
  NLITTsomf   = dNLITT * FLITTSOMF
  NLITTnmin   = dNLITT - NLITTsomf

  ! N SOM fast
  rNSOMF      = ((NSOMF*Runoff) / ROOTD) * RRUNBULK * 0.001
  dNSOMF      =  (NSOMF*dCSOMF) / CSOMF
  NSOMFsoms   = dNSOMF * FSOMFSOMS
  NSOMFnmin   = dNSOMF - NSOMFsoms

  ! N SOM slow
  dNSOMS      = (NSOMS*dCSOMS) / CSOMS

  ! N mineralisation, fixation, leaching, emission
  Nmineralisation = NLITTnmin + NSOMFnmin + dNSOMS
  Nfixation       = gCR * KNFIX 
  Nleaching       = (NMIN*RNLEACH*Drain) / WA
  Nemission       =  NMIN * KNEMIT * RWA
  fN2O            = 1. / (1. + exp(-RFN2O*(WFPS-WFPS50N2O)))
  NemissionN2O    = Nemission *     fN2O
  NemissionNO     = Nemission * (1.-fN2O) 

  end Subroutine CNsoil

  Subroutine Tsoil_calc
  integer iday,i
  Tsoil=0.
  do i=1,10
    Tsoil=Tsoil+Thist(i)/10.
  enddo
  do i=10,2,-1
    Thist(i)=Thist(i-1)
  enddo
!  Thist(1)=T(iday) 
  Thist(1)=T
  fTsoil = exp((Tsoil-10.)*(2.*TMAXF-Tsoil-10.)/(2.*TSIGMAF**2.)) 
  end Subroutine Tsoil_calc 
end module soil      
