Subroutine set_params(pa)

use parameters
implicit none
! As long as the total number of parameters stays below 100, the next line need not be changed
real pa(100)

! TREES
BETA      = pa(1)    ! (-)                        Sensitivity of LUE to [CO2]
CBTREE0   = pa(2)    ! (kg C tree-1)              Initial C in branches
CLTREE0   = pa(3)    ! (kg C tree-1)              Initial C in leaves
CRTREE0   = pa(4)    ! (kg C tree-1)              Initial C in roots
CSTREE0   = pa(5)    ! (kg C tree-1)              Initial C in tree stems
FB        = pa(6)    ! (-)                        Fraction of C allocated to branches
FLMAX     = pa(7)    ! (-)                        Max. fraction of C allocated to leaves
FNCLMIN   = pa(8)    ! (-)                        Min. N/C ratio of leaves as a fraction of the max.
FS        = pa(9)    ! (-)                        Fraction of C allocated to stems
FTCCLMIN  = pa(10)   ! (-)                        Min. longevity of leaves as a fraction of the max.
GAMMA     = pa(11)   ! (-)                        Autotrophic respiration as a fraction of GPP
KAC       = pa(12)   ! (m2)                       Allometric proportionality constant: crown area when branch C = 1 kg tree-1
KACEXP    = pa(13)   ! (-)                        Allometric exponent for crown area as a function of branch C
KBA       = pa(14)   ! (-)                        Proportionality constant for basal area as a function of tree volume/height ratio
KEXT      = pa(15)   ! (m2 m-2 leaf)              Light extinction coefficient
KH        = pa(16)   ! (m)                        Allometric proportionality constant: height when stem C = 1 kg tree-1
KHEXP     = pa(17)   ! (-)                        Allometric exponent for height as a function of stem C
KNMIN     = pa(18)   ! (kg N m-2)                 Michaelis-Menten Km-constant for uptake of mineral N
KNUPT     = pa(19)   ! (kg N kg-1 C d-1)          Michaelis-Menten Vmax-constant for uptake of mineral N
KRAININT  = pa(20)   ! ((mm d-1) (m2 leaf m-2)-1) Max. amount of rain interception per day as a function of LAI
LUEMAX    = pa(21)   ! (kg C MJ-1 PAR)            Max. Light-Use Efficiency at [CO2]=CO20
NCLMAX    = pa(22)   ! (kg N kg-1 C)              Max. N/C ratio of leaves
NCR       = pa(23)   ! (kg N kg-1 C)              N/C ratio of roots
NCW       = pa(24)   ! (kg N kg-1 C)              N/C ratio of wood
SLA       = pa(25)   ! (m2 leaf kg-1 C)           Specific Leaf Area
TCCB      = pa(26)   ! (d)                        Longevity of branches
TCCLMAX   = pa(27)   ! (d)                        Max. longevity of leaves
TCCR      = pa(28)   ! (d)                        Longevity of roots
TOPT      = pa(29)   ! (degC)                     Optimum temperature for LUE
TTOL      = pa(30)   ! (degC)                     Tolerance of LUE for suboptimal temperature
TRANCO    = pa(31)   ! (mm d-1)                   Transpiration coefficient (drought tolerance)
WOODDENS  = pa(32)   ! (kg C m-3 wood)            Wood density

! SOIL
CLITT0    = pa(33)   ! (kg C m-2)             Initial C in litter
CSOM0     = pa(34)   ! (kg C m-2)             Initial C in OM
CNLITT0   = pa(35)   ! (kg C kg-1 N)          Initial C/N ratio of litter
CNSOMF0   = pa(36)   ! (kg C kg-1 N)          Initial C/N ratio of fast-decomposing OM
CNSOMS0   = pa(37)   ! (kg C kg-1 N)          Initial C/N ratio of slowly decomposing OM
FCSOMF0   = pa(38)   ! (-)                    Initial C in fast-decomposing OM as a fraction of total OM
FLITTSOMF = pa(39)   ! (-)                    Fraction of decomposing litter that becomes OM
FSOMFSOMS = pa(40)   ! (-)                    Fraction of decomposing 'fast' OM that becomes slowly decomposing OM
RNLEACH   = pa(41)   ! (-)                    Mineral N concentration of drainage water as a ratio of that in soil water
KNEMIT    = pa(42)   ! (d-1)                  Max. relative emission rate of soil mineral N
NMIN0     = pa(43)   ! (kg N m-2)             Initial mineral N
TCLITT    = pa(44)   ! (d)                    Residence time of litter
TCSOMF    = pa(45)   ! (d)                    Residence time of fast-decomposing OM
TCSOMS    = pa(46)   ! (d)                    Residence time of slowly decomposing OM
TMAXF     = pa(47)   ! (degC)                 Temperature at which soil decomposition (fTsoil) is max.
TSIGMAF   = pa(48)   ! (degC)                 Tolerance of soil decomposition for suboptimal temperature
RFN2O     = pa(49)   ! (-)                    Sensitivity of the N2O/NO emission ratio to extreme values of water-filled pore space
WFPS50N2O = pa(50)   ! (-)                    Water-filled pore space at which the N2O and NO emission rates are equal

ROOTD     = pa(51)
WCST      = pa(52)
FWCAD     = pa(53)
FWCWP     = pa(54)
FWCFC     = pa(55)
FWCWET    = pa(56)

TREEDENS0 = pa(57)
LAT       = pa(58)

LAIMAX    = pa(59)  ! (m2 leaf m-2)

end Subroutine set_params
