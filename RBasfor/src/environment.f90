module environment

use parameters
implicit none
real    :: GR, PAR, T, RAIN, WN, VP, CO2A
integer :: doyi(NMAXDAYS), yeari(NMAXDAYS)
real    :: GRI(NMAXDAYS), TI(NMAXDAYS), RAINI(NMAXDAYS), WNI(NMAXDAYS), VPI(NMAXDAYS)
real    :: DAYL
real    :: Ndep

Contains

Subroutine set_weather_day(day, year,doy)
  integer :: day, doy, year
  integer :: y 
  doy   = doyi(day)
  GR    = GRI(day)
  PAR   = GR * 0.5
  T     = TI(day)
  RAIN  = RAINI(day)
  WN    = WNI(day)
  VP    = VPI(day)
  year  = yeari(day)
  y     = year - 1900
  CO2A  = 291 + 1.15*y - 0.03*(y**2) + 0.000354*(y**3) - 0.0000009*(y**4)
end Subroutine set_weather_day

Subroutine DDAYL(doy)
!=============================================================================
! Calculate day length (fractional days; -) from Julian day and latitude
!=============================================================================
  integer :: doy
  real    :: DEC, DECC, RAD
  RAD  = PI / 180.                                                    ! (radians deg-1)
  DEC  = -asin (sin (23.45*RAD)*cos (2.*PI*(doy+10.)/365.))           ! (radians)
  DECC = max(atan(-1./tan(RAD*LAT)),min( atan( 1./tan(RAD*LAT)),DEC)) ! (radians)
  DAYL = 0.5 * ( 1. + 2. * asin(tan(RAD*LAT)*tan(DECC)) / PI )        ! (d d-1)
end Subroutine DDAYL

Subroutine N_dep(year,doy,DAYS_NDEP,NDEPV)
  integer                  :: year,doy,j
  integer,dimension(100,2) :: DAYS_NDEP
  real   ,dimension(100  ) :: NDEPV
  integer                  :: idep
  real                     :: NDEPV_interval,t
  real   ,dimension(100)   :: tNdep
  t     = year           + (doy           -0.5)/366
  tNdep = DAYS_NDEP(:,1) + (DAYS_NDEP(:,2)-0.5)/366
  do j = 2,100
   if ( (tNdep(j-1)<t) .and. (tNdep(j)>=t) ) idep = j-1
  end do
  NDEPV_interval = NDEPV(idep+1) - NDEPV(idep)
  Ndep           = NDEPV(idep) + NDEPV_interval * (t            -tNdep(idep)) / &
                                                  (tNdep(idep+1)-tNdep(idep))
end Subroutine N_dep 

end module environment





