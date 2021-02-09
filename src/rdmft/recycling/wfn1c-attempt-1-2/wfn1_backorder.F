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
#include "global.fh"
#include "wfn1_param.fh"
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
      integer :: jj
      integer :: ilo, ihi, jlo, jhi
      integer :: iproc
      double precision, allocatable :: val(:)
c
      iproc = ga_nodeid()
      call ga_distribution(g_in,iproc,ilo,ihi,jlo,jhi)
      allocate(val(ilo:ihi))
      do jj = jlo, jhi
        call ga_get(g_in,ilo,ihi,jj,jj,val,ihi-ilo+1)
        call ga_put(g_out,ilo,ihi,icol(jj),icol(jj),val,ihi-ilo+1)
      enddo
      deallocate(val)
      end subroutine wfn1_backorder
c
c-----------------------------------------------------------------------
