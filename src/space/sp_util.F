!> \brief Find the lenght of the contents of a string
!>
!> @param[in] text  the input string
!> @return the length of the contents of "text"
!>
      integer function sp_strlen(text)
      implicit none
      character(len=*), intent(in) :: text
      integer :: ilen
      ilen = len(text)
      do while (ilen.gt.0.and.text(ilen:ilen).eq." ")
        ilen = ilen - 1
      enddo
      sp_strlen = ilen
      return
      end
!>
!> \brief Print IPL array
!>
!> IPL is a key data structure that lists information about the domain
!> composition. It contains a variety of data. Basically the array
!> consists of 3 sections:
!>   1. ipl(1:3,*) various sizes
!>   2. ipl(4:30,*) data needed to exchange atoms with neighboring boxes
!>                  (3x3x3 = 27 boxes including self)
!>   3. ipl(31:mbox,*) data about boxes that this processor owns
!> The data about locally owned boxes includes:
!>   - the box number
!>   - the range of solvent molecules in the box
!>   - the range of solute atoms in the box
!>   - etc.
!> This routine prints the contents of this array.
!>
!> @param[in] location  indicates the call location
!> @param[in] lfn       the file unit number
!>
      subroutine sp_pripl(location,lfn)
      implicit none
!
#include "mafdecls.fh"
#include "sp_common.fh"
!
      target :: int_mb
      integer, pointer :: ipl(:,:)
      character(len=*), intent(in)  :: location
      integer, intent(in) :: lfn
      integer, external :: sp_strlen
      integer :: ic
!
      ipl(1:mbox,1:mip2) => int_mb(i_ipl:i_ipl+mbox*mip2-1)
      write(lfn,1000)location(1:sp_strlen(location))
 1000 format('ipl from: ',a)
      write(lfn,1010)me,ipl(1,1),ipl(1,2),ipl(2,2)
 1010 format('proc: ',i8,/,'# boxes: ',i8,/,'# waters: ',i8,/,
     +       '# solute: ',i8)
      write(lfn,1020)"idx","boxnr","iwfr","iwto","isfr","isto"
 1020 format(a8,"  ",a8,"  ",2a8,"  ",2a8)
      do ic = 30+1, 30+ipl(1,1)
        write(lfn,1030)ic,ipl(ic,1),ipl(ic,2),ipl(ic,3),
     +                              ipl(ic,4),ipl(ic,5)
 1030   format(i8,"  ",i8,"  ",2i8,"  ",2i8)
      enddo
      write(lfn,1020)"num","jproc","iwfr","iwto","isfr","isto"
      do ic = 1, 27
        if (neighb(ic,2).ne.me) then
          write(lfn,1030)ic,neighb(ic,2),ipl(3+ic,1),ipl(3+ic,2),
     +                                   ipl(3+ic,3),ipl(3+ic,4)
        endif
      enddo
      end
