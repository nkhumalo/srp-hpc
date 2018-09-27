c-----------------------------------------------------------------------
c
      subroutine md_veladj(numw,nums,vw,vs,iwdx,isdx)
      implicit none
#include "md_common.fh"
#include "global.fh"
c
c     Adjust the velocities as specified in the input.
c
      integer numw ! number of waters on this processor
      integer nums ! number of solute atoms on this processor
      real*8  vw(mwm,3,mwa) ! water velocities
      real*8  vs(msa,3)     ! solute velocities
      integer iwdx(mwm)     ! local water molecule index
      integer isdx(msa,mis2)     ! local solute atom index
c
      integer isa
      integer ix
c
cDEBUG
c     write(*,*)'HVD: md_veladj'
cDEBUG
      if(ivelatm(1).le.0.or.ivelatm(2).le.0) return
cDEBUG
c     write(*,*)'HVD 2: md_veladj',ga_nodeid(),nums,ivelatm,mis2
cDEBUG
      do isa = 1, nums
cDEBUG
c     write(*,'("HVD 3: md_veladj",17i5)')ga_nodeid(),isa,
c    +      isdx(isa,1:mis2)
cDEBUG
      if(isdx(isa,1).ge.ivelatm(1).and.
     +   isdx(isa,1).le.ivelatm(2)) then
cDEBUG
c     write(*,*)'HVD 4: md_veladj',ivelatm,isa,isdx(isa,1)
cDEBUG
      do ix = 1,3
      vs(isa,ix)=vs(isa,ix)+veloadj(ix)
      enddo
      endif
      enddo
      ivelatm(1) = -ivelatm(1)
      ivelatm(2) = -ivelatm(2)
c
      end subroutine md_veladj
c
c-----------------------------------------------------------------------
