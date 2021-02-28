c-----------------------------------------------------------------------
c
      subroutine wfn1_diag(wfn1_param,alpha,g_x,g_t,g_exp)
      implicit none
C>
C> \brief Calculate diag(X_ii+alpha*X_ij) where X is a symmetric matrix
C>
C> NOTE: This routine is serial and does not synchronize. The rational
C>       for the latter is that in the application it is always followed
C>       by a ga_dgemm (which does synchronize). In addition where we
C>       invoke this is invoked 4 times where we use it, so where a
C>       single synchronization is sufficient to guarantee completion of
C>       all 4 matrices.
C>
C> Diagonalize a symmetric matrix where the off-diagonal elements have
C> been scaled with a factor alpha. The resulting eigenvector provide
C> a unitary transformation of the basis the matrix was expressed in.
C> Eigenvectors are sorted so the maximum coefficients appear on the 
C> diagonal to implement a maximum overlap ordering of the states.
C> 
#include "wfn1_param.fh"
#include "errquit.fh"
#include "global.fh"
!     Calculation parameters
      type(wfn1_prm), intent(in)   :: wfn1_param
c
      double precision, intent(in) :: alpha
c
      integer, intent(in)    :: g_x
      integer, intent(inout) :: g_t
      integer, intent(inout) :: g_exp
c
      double precision :: dnormX  ! The norm of alpha*X
      double precision :: dnormXn ! The norm of the n-th term in the
                                  ! Taylor series of exp(alpha*X)
c
      double complex, allocatable :: skew(:,:) ! skew symmetric mat
c     double complex, allocatable :: skvl(:,:) ! left eigenvectors
      double complex              :: skvl      ! left eigenvectors
      double complex, allocatable :: skvr(:,:) ! right eigenvectors
      double complex, allocatable :: eig(:)    ! eigenvalues
c
      double precision, allocatable :: rwork(:)   ! scratch space
c
c     double precision, allocatable, target :: buf(:)
      double precision, allocatable :: block(:,:)
c
      double precision, allocatable :: u(:,:) ! the unitary matrix
      double precision              :: dnrm
      double precision, external    :: ddot
c
      integer, allocatable :: icol(:)
c
      integer :: itypx, nrowx, ncolx
      integer :: info ! status flag returned by LAPACK
c
      integer :: iproc, nproc
      integer :: ip, ii, jj, kk ! counters
      integer :: ixlo, ixhi, jxlo, jxhi ! local block limits
      integer :: itlo, ithi, jtlo, jthi ! local block limits
      integer :: ixdim, jxdim
      integer :: itdim, jtdim
c
      logical :: o_repeat
c
      character*9 pname
      parameter(pname="wfn1_diag:")
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call ga_inquire(g_x,itypx,nrowx,ncolx)
      if (nrowx.ne.ncolx) then
        call errquit(pname//" matrix g_x not square",0,UERR)
      endif
      call ga_distribution(g_x,iproc,ixlo,ixhi,jxlo,jxhi)
      ixdim = ixhi-ixlo+1
      jxdim = jxhi-jxlo+1
      call ga_distribution(g_t,iproc,itlo,ithi,jtlo,jthi)
      itdim = ithi-itlo+1
      jtdim = jthi-jtlo+1
c
      allocate(block(ixlo:ixhi,jxlo:jxhi))
      allocate(eig(1:nrowx))
      allocate(icol(1:nrowx))
c
      call ga_get(g_x,ixlo,ixhi,jxlo,jxhi,block,ixdim)
      do jj = jxlo, jxhi
        do ii = ixlo, ixhi
          if (ii.ne.jj) then
            block(ii,jj) = alpha*block(ii,jj)
          endif
        enddo
      enddo
      call ga_put(g_exp,ixlo,ixhi,jxlo,jxhi,block,ixdim)
      call ga_sync()
c
      call ga_diag_std(g_exp,g_t,eig)
      call wfn1_column_max(wfn1_param,g_t,icol)
      call wfn1_backorder(wfn1_param,g_t,icol,g_exp)
c
      deallocate(icol)
      deallocate(block)
      deallocate(eig)
c
      end
c
c-----------------------------------------------------------------------
