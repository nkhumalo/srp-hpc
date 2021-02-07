c-----------------------------------------------------------------------
c
      subroutine wfn1_diag(wfn1_param,alpha,g_x,g_t,g_exp,eig)
      implicit none
C>
C> \brief Calculate diag(X_ii+alpha*X_ij) where X is a hermitian matrix
C>
C> Diagonalize a hermitian matrix where the off-diagonal elements have
C> been scaled with a factor alpha. The resulting eigenvector provide
C> a unitary transformation of the basis the matrix was expressed in.
C> Eigenvectors are sorted so the maximum coefficients appear on the 
C> diagonal to implement a maximum overlap ordering of the states.
C> 
!include "wfn1_param.fh"
#include "errquit.fh"
#include "global.fh"
!     Calculation parameters
      type(wfn1_prm), intent(in)   :: wfn1_param
c
      double precision, intent(in) :: alpha
c
!>    the Hermitian matrix to diagonalize
      integer, intent(in) :: g_x
!>    temporary
      integer, intent(in) :: g_t
!>    the resulting eigen vectors sorted by maximum overlap
      integer, intent(in) :: g_exp
!>    the eigen values sorted by maximum overlap
      double precision, intent(out) :: eig(1:wfn1_param%nmo)
c
      double precision :: dnormX  ! The norm of alpha*X
      double precision :: dnormXn ! The norm of the n-th term in the
                                  ! Taylor series of exp(alpha*X)
c
      double complex, allocatable :: skew(:,:) ! skew symmetric mat
c     double complex, allocatable :: skvl(:,:) ! left eigenvectors
      double complex              :: skvl      ! left eigenvectors
      double complex, allocatable :: skvr(:,:) ! right eigenvectors
      double complex, allocatable :: deig(:)   ! eigenvalues
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
      integer :: itypt, nrowt, ncolt
      integer :: itype, nrowe, ncole
      integer :: info ! status flag returned by LAPACK
c
      integer :: iproc, nproc
      integer :: ip, ii, jj, kk ! counters
      integer :: ixlo, ixhi, jxlo, jxhi ! local block limits
      integer :: itlo, ithi, jtlo, jthi ! local block limits
      integer :: ixdim, jxdim
      integer :: itdim, jtdim
      integer :: nb
c
      logical :: o_repeat
c
      character*9 pname
      parameter(pname="wfn1_diag:")
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
#ifndef SCALAPACK
You must have ScaLAPACK!
#endif
      call ga_inquire(g_x,itypx,nrowx,ncolx)
      if (nrowx.ne.ncolx) then
        call errquit(pname//" matrix g_x not square",0,UERR)
      endif
      call ga_inquire(g_t,itypt,nrowt,ncolt)
      call ga_inquire(g_exp,itype,nrowe,ncole)
      if (itypx.ne.itypt.or.itypx.ne.itype) then
        call errquit(pname//" matrices types mismatch",0,UERR)
      endif
      allocate(deig(1:wfn1_param%nmo))
      allocate(icol(1:wfn1_param%nmo))
      call ga_copy(g_x,g_exp)
      call wfn1_scale_offdiag(g_exp,alpha)
      nb = 0
      if (itypx.eq.MT_DBL) then
        call ga_pdsyevd(g_exp,g_t,deig,nb)
      else if (itypx.eq.MT_DCPL) then
        call ga_pzheevd(g_exp,g_t,deig,nb)
      else
        call errquit(pname//" invalid type",itypx,UERR)
      endif
      call ga_diag_std(g_exp,g_t,deig)
      call wfn1_column_max(wfn1_param,g_t,icol)
      call wfn1_backorder(wfn1_param,g_t,icol,g_exp)
      call wfn1_backorder_e(wfn1_param,deig,icol,eig)
      deallocate(icol)
      deallocate(deig)
c
      end
c
c-----------------------------------------------------------------------
