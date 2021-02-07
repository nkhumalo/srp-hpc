c-----------------------------------------------------------------------
c
C> \brief Order the vectors back
C>
C> The diagonalization orders the eigenvectors in increasing order of
C> eigenvalues. This routine moves them back to their original ordering.
C>
C> The array `icol` lists where every vector came from, we will use this
C> information to put the vector back to where it came from. 
C>
      subroutine wfn1_backorder(wfn1_param,g_in,icol,g_out)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!include "wfn1_param.fh"
!     The calculation parameters
      type(wfn1_prm), intent(in) :: wfn1_param
!     The input vectors
      integer, intent(in) :: g_in
!     The columns where the vectors came from
      integer, intent(in) :: icol(1:wfn1_param%nmo)
!     The output vectors
      integer, intent(inout) :: g_out
c
c     Local
c
      integer :: ii,jj,itypei,itypeo
      integer :: ilo, ihi, jlo, jhi
      integer :: iproc
      double precision, allocatable :: val(:)
      complex(kind=8),  allocatable :: zval(:)
c
      iproc = ga_nodeid()
      call ga_inquire(g_in,itypei,ii,jj)
      call ga_inquire(g_out,itypeo,ii,jj)
      if (itypei.ne.itypeo) then
        call errquit("wfn1_backorder: type mismatch",0,UERR)
      endif
      call ga_distribution(g_in,iproc,ilo,ihi,jlo,jhi)
      if (itypei.eq.MT_DBL) then
        allocate(val(ilo:ihi))
        do jj = jlo, jhi
          call ga_get(g_in,ilo,ihi,jj,jj,val,ihi-ilo+1)
          call ga_put(g_out,ilo,ihi,icol(jj),icol(jj),val,ihi-ilo+1)
        enddo
        deallocate(val)
      else if (itypei.eq.MT_DCPL) then
        allocate(zval(ilo:ihi))
        do jj = jlo, jhi
          call ga_get(g_in,ilo,ihi,jj,jj,zval,ihi-ilo+1)
          call ga_put(g_out,ilo,ihi,icol(jj),icol(jj),zval,ihi-ilo+1)
        enddo
        deallocate(zval)
      else
        call errquit("wfn1_backorder: invalid type",itypei,UERR)
      endif
      end subroutine wfn1_backorder
c
c-----------------------------------------------------------------------
c
C> \brief Order the eigenvalues back
C>
C> The diagonalization orders the eigenvalues in increasing order of
C> eigenvalues. This routine moves them back to their original ordering.
C>
C> The array `icol` lists where every eigenvalue came from, we will use
C> this information to put the vector back to where it came from. 
C>
      subroutine wfn1_backorder_e(wfn1_param,e_in,icol,e_out)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!include "wfn1_param.fh"
!     The calculation parameters
      type(wfn1_prm), intent(in) :: wfn1_param
!     The input vectors
      double precision, intent(in) :: e_in(1:wfn1_param%nmo)
!     The columns where the vectors came from
      integer, intent(in) :: icol(1:wfn1_param%nmo)
!     The output vectors
      double precision, intent(out) :: e_out(1:wfn1_param%nmo)
c
c     Local
c
      integer :: ii,jj,itypei,itypeo
      integer :: ilo, ihi, jlo, jhi
      integer :: iproc
      double precision, allocatable :: val(:)
      complex(kind=8),  allocatable :: zval(:)
c
      e_out = +1.0d300
      do jj = 1, wfn1_param%nmo
        e_out(icol(jj)) = e_in(jj)
      enddo
      end subroutine wfn1_backorder_e
c
c-----------------------------------------------------------------------
