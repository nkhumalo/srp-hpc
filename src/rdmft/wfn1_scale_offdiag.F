c-----------------------------------------------------------------------
c
C> \brief Scale the off-diagonal elements of a matrix by the given
C> factor
C>
C> When using diagonalization based methods it is sometimes necessary
C> to reduce the rotations the diagonalization imposes. One way to do
C> this is to reduce the off-diagonal matrix elements by some factor.
C> This routine takes a matrix, scales the off-diagonal elements and
C> returns the resulting matrix.
C>
      subroutine wfn1_scale_offdiag(g_a,fac)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!     The matrix 
      integer,       intent(inout) :: g_a
!     The scale factor for the off-diagonal elements
      double precision, intent(in) :: fac
c
c     Local
c
      integer :: iproc
      integer :: ilo, ihi
      integer :: jlo, jhi
      integer :: ii,  jj
      integer :: itype
c
      double precision, allocatable :: val(:,:)
      complex(kind=8),  allocatable :: zval(:,:)
c
      iproc = ga_nodeid()
      call ga_inquire(g_a,itype,ii,jj)
      call ga_distribution(g_a,iproc,ilo,ihi,jlo,jhi)
      if (itype.eq.MT_DBL) then
        allocate(val(ilo:ihi,jlo:jhi))
        call ga_get(g_a,ilo,ihi,jlo,jhi,val,ihi-ilo+1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            if (ii.ne.jj) then
              val(ii,jj) = fac*val(ii,jj)
            endif
          enddo
        enddo
        call ga_put(g_a,ilo,ihi,jlo,jhi,val,ihi-ilo+1)
        deallocate(val)
      else if (itype.eq.MT_DCPL) then
        allocate(zval(ilo:ihi,jlo:jhi))
        call ga_get(g_a,ilo,ihi,jlo,jhi,zval,ihi-ilo+1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            if (ii.ne.jj) then
              zval(ii,jj) = fac*zval(ii,jj)
            endif
          enddo
        enddo
        call ga_put(g_a,ilo,ihi,jlo,jhi,zval,ihi-ilo+1)
        deallocate(zval)
      else
        call errquit("wfn1_scale_offdiag: invalid type",itype,UERR)
      endif
      call ga_sync()
c
      end subroutine wfn1_scale_offdiag
c
c-----------------------------------------------------------------------
