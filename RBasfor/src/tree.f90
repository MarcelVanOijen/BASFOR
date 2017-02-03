Module tree

use parameters
use management
use environment
implicit none

! Seasonality
integer :: treegrow, leaffall
real    :: dchillday, dtsum
! Morphology
real :: AC, ACtree, BA, BAtree, CBtree, CStree, DBH, H, LA, LAI
! foliarDynamics
real :: dCL, dCLold, dCLprun, dCLsen, dCLsenNdef, dCLthin, dCRESgrow, dCRESthin
real :: dLAIold, LAIsurv
real :: dNL, dNLdeath, dNLlitt, dNLthin
real :: NLMax, NLMin, NLsurvMax, NLsurvMin
real :: NEFF, NLsurv, NSTATUS, NSTATUSF, recycNLold, retrNLmax
! NPP
real :: fLUECO2, fLUET, GPP, LUE, NPPmaxN, PARabsCrown, PARabs
! Allocation
real :: FL, FR
! NdemandOrgans 
real :: gCBmaxN , gCLmaxN , gCRmaxN , gCSmaxN, gCRESmaxN
real :: NdemandB, NdemandL, NdemandR, NdemandS, NdemandRES
! gtreeNupt
real :: fNgrowth, gCB, gCL, gCR, gCS, gCRES, gNB, gNL, gNR, gNS, gNRES
real :: Ndemand, Nupt
real :: retrNL
! CNtree
real :: dCB, dCBprun, dCBsen, dCBthin
real :: dCR,          dCRsen, dCRthin
real :: dCS
real :: dNBlitt, dNRsomf
real :: harvCS, harvNS, LAIcrit

Contains

  Subroutine dtsum_dchillday(chillday,doy,Tsum, dchillday,dTsum)
    integer            :: doy
    real               :: chillday , Tsum
    real               :: dchillday, dTsum
    integer, parameter :: chillstart = 304, Tsumstart = 32
#ifdef deciduous
!   DECIDUOUS TREES
    dchillday = 0
    dTsum     = 0
    if (doy==chillstart) then
      dchillday = -chillday
    else if (T<Tc) then
      dchillday = 1
    end if
    if (doy<Tsumstart) then 
      dTsum = -Tsum
    else if (T>Tb) then
      dTsum = T-Tb
    end if
#endif
  end Subroutine dtsum_dchillday  

  Subroutine phenology(chillday,doy,Tsum, leaffall,treegrow)
    integer :: doy 
    real    :: chillday, Tsum
    integer :: leaffall, treegrow
#ifdef deciduous
!   DECIDUOUS TREES
    integer :: iautumn
    real    :: Tcrit
!   Winter
    iautumn  = 0
    leaffall = 0 
    treegrow = 0
!   Autumn
    call DDAYL(doy)
    if ( (doy>200) .and. (DAYL<DAYLAUT) ) then
      iautumn  = 1
      leaffall = 1 
    end if  
!   Spring  
    Tcrit = Tcrita - Tcritb * max(0.0,LOG(chillday))
    if ( (Tsum>Tcrit) .and. (iautumn/=1) ) treegrow = 1
#else
!   CONIFEROUS TREES
    leaffall = 1
    treegrow = 1
#endif
  end Subroutine phenology  

  Subroutine morphology(CL,CS,CB,treedens, LAI)
    real :: CB, CL, CS, LAI, treedens
    CBtree    = CB / treedens
    CStree    = CS / treedens
    H         = KH  * (CStree**KHEXP)                              ! m
    BAtree    = KBA * (CStree/WOODDENS) / H                        ! m2 tree-1
    BA        = BAtree * treedens                                  ! m2 m-2
    DBH       = 2 * (BAtree/PI)**0.5                               ! m
    ACtree    = KAC * (CBtree**KACEXP)                             ! m2 AC tree-1
    AC        = ACtree * treedens
      if (AC.gt.1.0) AC = 1.0                                      ! m2 AC m-2
    LA        = CL * SLA                                           ! m2 leaf m-2
    LAI       = LA / AC                                            ! m2 leaf m-2 AC
  end Subroutine morphology  

  Subroutine foliarDynamics(CL,CRES,fTran,NL,LAI)
    real :: CL, CRES, LAI, fTran, NL
  ! Leaf N
    NLMAX     = NLAMAX *                LAI                        ! kg N m-2 AC
    NLMIN     = NLAMAX * (1.0-exp(-KEXT*LAI)) / KEXT               ! kg N m-2 AC
    if (LAI == 0.) then
      NEFF    = 0.
      NSTATUS = 0.
    else
      NEFF    = max( 0., min( 1., (NL/AC-NLAMIN*LAI) / (NLMAX-NLAMIN*LAI) ) ) ! -
      NSTATUS = max( 0., (NL/AC-NLMIN)/(NLMAX-NLMIN) )             ! -
    end if
    NSTATUSF  = min( NSTATUS, 0.999*FNCLMIN )                      ! -
    LAIcrit   = log( (NSTATUSF-1.0)/(NSTATUSF-FNCLMIN) ) / KEXT    ! m2 leaf m-2 AC  
  ! Leaf C
    dCLthin   = CL   * thinFRT                                     ! kg C m-2 d-1
    dCRESthin = CRES * thinFRT
    dCLprun   = CL   * prunFRT                                     ! kg C m-2 d-1
  ! Reserves are assumed to be released when we are in the growing season
    dCLsen    = leaffall * ( CL / (fTran+FTCCLMIN*(1.-fTran)) ) / TCCLMAX
    dCRESgrow = treegrow * CRES / DELT
    dNLthin   = NL * thinFRT                                       ! kg N m-2 d-1
    if (thinFRT==0) then
      dCLsenNdef = ( (LAI-LAIcrit)/DELT ) * AC / SLA               ! kg C m-2 d-1
      dCLold     = max( dCLsenNdef, dCLprun+dCLsen )               ! kg C m-2 d-1
      dLAIold    = dCLold * SLA / AC                               ! m2 leaf m-2 AC d-1
      LAIsurv    = LAI - dLAIold * DELT                            ! m2 leaf m-2 AC
      NLsurvMAX  = NLAMAX *                LAIsurv                 ! kg N m-2 AC
      NLsurvMIN  = NLAMAX * (1.0-exp(-KEXT*LAIsurv)) / KEXT        ! kg N m-2 AC
      NLsurv     = NLsurvMIN + (NLsurvMAX-NLsurvMIN)*NSTATUS       ! kg N m-2 AC
      dNLdeath   = (NL-NLsurv*AC) / DELT                           ! kg N m-2 d-1
      dNLlitt    = dLAIold * NLAMIN * AC                           ! kg N m-2 d-1
      recycNLold = dNLdeath - dNLlitt                              ! kg N m-2 d-1
      retrNLMAX  = (NLsurv-NLsurvMIN) * AC / DELT                  ! kg N m-2 d-1
    else
      dCLsen     = 0
      dCLsenNdef = 0
      dCLold     = 0
      dLAIold    = 0
      LAIsurv    = LAI
      NLsurvMAX  = NL
      NLsurvMIN  = NL
      NLsurv     = NL
      dNLdeath   = 0
      dNLlitt    = dNLthin
      recycNLold = 0
      retrNLMAX  = 0
    end if  
    dCL       = dCLthin + dCLold                                     ! kg C m-2 d-1
    dNL       = dNLthin + dNLdeath                                   ! kg N m-2 d-1

    end Subroutine foliarDynamics  

  Subroutine NPP(fTran)
    integer :: day
    real    :: fTran
    fLUECO2     = 1.0 + BETA * log(CO2A/CO20)
    fLUET       = exp( -0.5*((T - TOPT)/TTOL)**2. )
    LUE         = LUEMAX * fLUECO2 * fLUET * fTran * NEFF
    PARabsCrown = PAR * (1. - exp(-KEXT*LAI))
    PARabs      = PARabsCrown * AC
    GPP         = PARabs * LUE
    NPPmaxN     = GPP * (1.-GAMMA)
  end Subroutine NPP

  Subroutine allocation(fTran,LAI)
  real :: fTran,LAI
    FL = FLMAX * fTran * NEFF                                     ! -
    if (LAI>LAIMAX) FL = 0
    FR = (1 - (FL+FB+FS))
  end Subroutine allocation

  Subroutine NdemandOrgans
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2011-10-05
!!! The sinks for the growing organs (L,B,S,R) are assumed to be active only
!!! during the growing season, whereas the reserves-sink takes over during fall.
!!! The nitrogen demanded for reserves is proportional to NPP, with an N-C ratio
!!! equal to NCLMAX, so reserves are ready for growth of young, high-N leaves.
    gCLmaxN    = FL * NPPmaxN *    treegrow
    gCBmaxN    = FB * NPPmaxN *    treegrow
    gCSmaxN    = FS * NPPmaxN *    treegrow
    gCRmaxN    = FR * NPPmaxN *    treegrow
    gCRESmaxN  =      NPPmaxN * (1-treegrow)
  ! 2010-08-10: New foliar dynamics
    NdemandL   = gCLmaxN   * NCLMAX                               ! kg N m-2 d-1
    NdemandB   = gCBmaxN   * NCW                                  ! kg N m-2 d-1
    NdemandS   = gCSmaxN   * NCW                                  ! kg N m-2 d-1
    NdemandR   = gCRmaxN   * NCR                                  ! kg N m-2 d-1
    NdemandRES = gCRESmaxN * NCLMAX                               ! kg N m-2 d-1
  end Subroutine NdemandOrgans

  Subroutine gtreeNupt(Nsup)
    real Nsup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 2011-10-05
!!! Soil nitrogen is assumed to be available to the plants only during the
!!! growing season, so thereafter soil N-supply becomes zero.
!!! During the fall, N-dermand of reserves is assumed to be proportional to NPP
!!! (which in that season only has the reserves as sink) whereas N-supply is
!!! from recycling of old leaves only (soil N no longer being available).
    real Nsupgrow  
    Nsupgrow = Nsup * treegrow                                    ! kg N m-2 d-1
    Ndemand  = NdemandL+NdemandB+NdemandS+NdemandR + NdemandRES   ! kg N m-2 d-1
    if  (Ndemand < (recycNLold + Nsupgrow)) then
      Nupt   = Ndemand - recycNLold                               ! kg N m-2 d-1
      retrNL = 0                                                  ! kg N m-2 d-1
    else if (Ndemand < recycNLold + Nsupgrow + retrNLMAX) then
      Nupt   = Nsupgrow                                           ! kg N m-2 d-1
      retrNL = Ndemand - recycNLold - Nsupgrow                    ! kg N m-2 d-1
    else
      Nupt   = Nsupgrow                                           ! kg N m-2 d-1
      retrNL = retrNLMAX                                          ! kg N m-2 d-1
    end if
    if (Ndemand > 0.) then
      fNgrowth = max( 0., min( 1., (recycNLold+Nsupgrow+retrNLMAX) / Ndemand ) )
    else
      fNgrowth = 0.
    end if                                                        ! -
    gNL      = NdemandL   * fNgrowth                              ! kg N m-2 d-1
    gNB      = NdemandB   * fNgrowth                              ! kg N m-2 d-1
    gNS      = NdemandS   * fNgrowth                              ! kg N m-2 d-1
    gNR      = NdemandR   * fNgrowth                              ! kg N m-2 d-1
    gNRES    = NdemandRES * fNgrowth                              ! kg N m-2 d-1
    gCL      = gNL   / NCLMAX                                     ! kg C m-2 d-1
    gCB      = gNB   / NCW                                        ! kg C m-2 d-1
    gCS      = gNS   / NCW                                        ! kg C m-2 d-1
    gCR      = gNR   / NCR                                        ! kg C m-2 d-1
    gCRES    = gNRES / NCLMAX                                     ! kg C m-2 d-1
  end Subroutine gtreeNupt

  Subroutine CNtree(CR,CS,CB)
    real :: CB, CR, CS
  ! Branch C
    dCBthin    = CB * thinFRT                                     ! kg C m-2 d-1
    dCBprun    = CB * prunFRT                                     ! kg C m-2 d-1
    dCBsen     = CB / TCCB                                        ! kg C m-2 d-1
    dCB        = dCBthin + dCBprun + dCBsen                       ! kg C m-2 d-1
    dNBlitt    = dCB * NCW                                        ! kg N m-2 d-1
  ! Stem C & harvest
    dCS        = CS * thinFRT                                     ! kg C m-2 d-1
    harvCS     = dCS + dCRESthin                                  ! kg C m-2 d-1
    harvNS     = dCS * NCW                                        ! kg N m-2 d-1
  ! Roots C
    dCRthin    = CR * thinFRT                                     ! kg C m-2 d-1
    dCRsen     = CR / TCCR                                        ! kg C m-2 d-1
    dCR        = dCRthin + dCRsen                                 ! kg C m-2 d-1
    dNRsomf    = dCR * NCR                                        ! kg N m-2 d-1
  end Subroutine CNtree
  
end Module tree 
