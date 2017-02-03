module parameters

! Run control and mathematical constants
  integer, parameter :: DELT  = 1
  real   , parameter :: PI    = ACOS(-1.0)
  
! Environment
  integer, parameter :: NMAXDAYS = 60000
  real   , parameter :: CO20  = 350    ! (ppm) Reference value of [CO2] at which fLUECO2 = 1
  real               :: LAT

! Soil conditions 
  real, parameter    :: KNFIX = 0, KRUNOFF = 0.5, RRUNBULK = 0.05, SLOPE = 0
  real               :: FWCAD, FWCWP, FWCFC, FWCWET, WCST
  real               :: ROOTD
 
! Forest management
  real               :: TREEDENS0

! TREES
real :: BETA, CBTREE0, CLTREE0, CRTREE0, CSTREE0, FB, FLMAX, FNCLMIN, FS
real :: FTCCLMIN, GAMMA, KAC, KACEXP, KBA, KEXT, KH, KHEXP, KNMIN, KNUPT
real :: KRAININT, LAI0, LUEMAX, NCLMAX, NCR, NCW, SLA, TCCB, TCCLMAX, TCCR, TOPT
real :: TTOL, TRANCO, WOODDENS
real :: LAIMAX

! Composite parameters
real :: NLAMAX, NLAMIN

! SOIL
real :: CLITT0, CSOM0, CNLITT0, CNSOMF0, CNSOMS0, FCSOMF0, FLITTSOMF, FSOMFSOMS
real :: RNLEACH, KNEMIT, NMIN0, TCLITT, TCSOMF, TCSOMS, TMAXF, TSIGMAF, RFN2O
real :: WFPS50N2O

! Phenological parameters determined by calibration using the average temperature (DCAM)
#ifdef deciduous
real, parameter :: DAYLAUT   = 0.4498876
real, parameter :: Tb        = 9.2426324
real, parameter :: Tc        = -3.3501598 ! (degC d)
real, parameter :: Tcrita    = 35.5237816
real, parameter :: Tcritb    = 1.7617968
#endif
  
end module parameters
