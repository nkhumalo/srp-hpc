c-----------------------------------------------------------------------
c
C> \brief Change the sign of the off-diagonal elements
C>
C> If the matrix is an orthogonal transformation (e.g. a rotation)
C> then this operation changes the direction of the transformation.
C>
      subroutine wfn1_sgn_offdiag(g_momo)
      implicit none
c
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
!include "wfn1_param.fh"
!include "wfn1_drv.fh"
c
!>    The matrix
      integer,        intent(inout) :: g_momo
c
!>    My processor rank
      integer :: iproc
!>    Lower and upper block limits
      integer :: ilo, ihi, jlo, jhi
!>    i-dimension, and j-dimension
      integer :: ni, nj
!>    Counters
      integer :: ii, jj
!>    Data type
      integer :: itype
c
!>    Buffers
      double precision, allocatable :: buf(:,:)
      complex(kind=8),  allocatable :: zbuf(:,:)
c
      iproc = ga_nodeid()
c
      call ga_sync()
      call ga_inquire(g_momo,itype,ni,nj)
      call ga_distribution(g_momo,iproc,ilo,ihi,jlo,jhi)
c
      if (itype.eq.MT_DBL) then
c
        allocate(buf(ilo:ihi,jlo:jhi))
        call ga_get(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            if (ii.ne.jj) then
              buf(ii,jj) = -buf(ii,jj)
            endif
          enddo ! ii
        enddo ! jj
        call ga_put(g_momo,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
        deallocate(buf)
c
      else if (itype.eq.MT_DCPL) then
c
        allocate(zbuf(ilo:ihi,jlo:jhi))
        call ga_get(g_momo,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            if (ii.ne.jj) then
              zbuf(ii,jj) = -zbuf(ii,jj)
            endif
          enddo ! ii
        enddo ! jj
        call ga_put(g_momo,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
        deallocate(zbuf)
c
      else
c
        call errquit("wfn1_sgn_offdiag: invalid data type",itype)
c
      endif
c
      call ga_sync()
c
      end subroutine wfn1_sgn_offdiag
