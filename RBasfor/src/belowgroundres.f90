module belowgroundres

use parameters
use management
use environment
implicit none

real :: BBRAD, NRADC, NRADS, RLWN
real :: RAINint, SSLOPE, SVP, WDF
real :: BOLTZM, LHVAP, PSYCH
real :: PENMD, PENMRC, PENMRS, PEvap, PTran

real :: fAvail, FRR, RWA, WAAD, WC, WCAD, WCCR, WCFC, WCWP, WCWET, WFPS
real :: Evap, Tran, fTran

real :: Nsup

Contains

  Subroutine PETtr(LAI,RAINint)
  !=============================================================================
  ! Calculate potential rates of evaporation and transpiration (mm d-1)
  ! Author - Marcel van Oijen (CEH-Edinburgh)
  !=============================================================================
  integer :: day
  real    :: LAI
  real    :: RAINint
  
  RAINint = min(RAIN,KRAININT*LAI)                         ! % (mm d-1)
      
  BOLTZM = 5.668E-8                                        ! % (J m-2 s-1 K-4)
  LHVAP  = 2.4E6                                           ! % (J kg-1)
  PSYCH  = 0.067                                           ! % (kPA degC-1))
  BBRAD  = BOLTZM * (T+273.)**4 * 86400.                   ! % (J m-2 d-1)
  SVP    = 0.611 * exp(17.4 * T/(T+239.))                  ! % (kPa)
  SSLOPE  = 4158.6 * SVP / (T+239.)**2                     ! % (kPA degC-1)
  RLWN   = BBRAD * max(0.,0.55*(1.-VP/SVP))                ! % (J m-2 d-1)
  NRADS  = GR*1.E6 * (1.-0.15) - RLWN                      ! % (J m-2 d-1)
  NRADC  = GR*1.E6 * (1.-0.25) - RLWN                      ! % (J m-2 d-1)
  PENMRS = NRADS * SSLOPE/(SSLOPE+PSYCH)                   ! % (J m-2 d-1)
  PENMRC = NRADC * SSLOPE/(SSLOPE+PSYCH)                   ! % (J m-2 d-1)
  WDF    = 2.63 * (1.0 + 0.54 * WN)                        ! % (kg m-2 d-1 kPa-1)
  PENMD  = LHVAP * WDF * (SVP-VP) * PSYCH/(SSLOPE+PSYCH)   ! % (J m-2 d-1)

  PEvap  =     exp(-0.5*LAI)  * (PENMRS + PENMD) / LHVAP   ! % (mm d-1)
  PTran  = (1.-exp(-0.5*LAI)) * (PENMRC + PENMD) / LHVAP   ! % (mm d-1)
  PTran  = max( 0., PTran-0.5*RAINint )                    ! % (mm d-1)

  end Subroutine PETtr

  Subroutine water_flux(WA,Evap,Tran,fTran,RWA,WFPS)
  !=============================================================================
  ! Calculate rates of evaporation and transpiration (mm d-1), and the
  ! transpiration realisation factor (-)
  ! Author: Marcel van Oijen (CEH-Edinburgh)
  ! Date:   6-11-2005, Modified: 10-7-2012
  !=============================================================================
  real :: WA
  real :: Evap, Tran, fTran, RWA, WFPS
  real :: RWAevap
  
  WCAD   = WCST * FWCAD                                          ! % (m3 m-3)
  WCWP   = WCST * FWCWP                                          ! % (m3 m-3)
  WCFC   = WCST * FWCFC                                          ! % (m3 m-3)
  WCWET  = WCST * FWCWET                                         ! % (m3 m-3)

  WC     = 0.001 * WA   / ROOTD                                  ! % (m3 m-3)
  RWAevap = max(0., min(1., (WC - WCWP) / (WCWET - WCWP) ) )     ! % (-)
  RWA     = max(0., min(1., (WC - WCAD) / (WCFC  - WCAD) ) )     ! % (-)
  WFPS    = max(0., min(1., (WC - WCAD) / (WCST  - WCAD) ) )     ! % (-)
  WAAD   = 1000. * WCAD * ROOTD                                  ! % (mm)
      
  Evap   = PEvap * RWAevap                                       ! % (mm d-1)
  WCCR = WCWP + max( 0.01, PTran/(PTran+TRANCO) * (WCFC-WCWP) )  ! % (m3 m-3)
  if (WC.gt.WCCR) then
       FRR = max(0., min(1., (WCST-WC)/(WCST-WCWET) ))
  else
       FRR = max(0., min(1., (WC-WCWP)/(WCCR-WCWP)  ))
  endif                                                          ! % (mm mm-1)
  Tran   = PTran * FRR                                           ! % (mm d-1)
  if (Evap+Tran.gt.0) then
       fAvail = min( 1., ((WA-WAAD)/DELT) / (Evap+Tran) )
  else
       fAvail = 0.                                                
  endif                                                          ! % (mm mm-1)
         
  Evap   = Evap * fAvail                                         ! % (mm d-1)
  Tran   = Tran * fAvail                                         ! % (mm d-1)

  if (PTran.gt.0) then
     fTran = Tran / PTran                                        ! % (-)
  else
     fTran = 1.                                                  ! % (-)
  endif
    
  end Subroutine water_flux
  
  Subroutine Nsupply(CR,NMIN,Nsup)
  real :: CR
  real :: NMIN
  real :: Nsup
  Nsup = min(CR*KNUPT*(NMIN/(KNMIN+NMIN)),NMIN/DELT)
  end Subroutine Nsupply
  
end module belowgroundres      
