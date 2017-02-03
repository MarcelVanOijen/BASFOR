module management

use parameters
use environment
implicit none
real    :: Nfert, prunFRT, thinFRT, thintreedens, treedens

Contains     

  Subroutine fert_prune_thin(year,doy,DAYS_FERT ,NFERTV , &
                                      DAYS_PRUNT,FRPRUNT, &
                                      DAYS_THINT,FRTHINT)
  integer                  :: year,doy,i
  integer,dimension(100,2) :: DAYS_FERT, DAYS_PRUNT, DAYS_THINT
  real   ,dimension(100  ) :: NFERTV   , FRPRUNT   , FRTHINT
  Nfert   = 0
  prunFRT = 0
  thinFRT = 0
  do i=1,100    
    if ( (year==DAYS_FERT (i,1)) .and. (doy==DAYS_FERT (i,2)) ) then
      Nfert   = NFERTV (i)
	end if
    if ( (year==DAYS_PRUNT(i,1)) .and. (doy==DAYS_PRUNT(i,2)) ) then
      prunFRT = FRPRUNT(i)
	end if
    if ( (year==DAYS_THINT(i,1)) .and. (doy==DAYS_THINT(i,2)) ) then
      thinFRT = FRTHINT(i)
	end if
  end do
  thintreedens = treedens * thinFRT
  end Subroutine fert_prune_thin
    
end module management      
