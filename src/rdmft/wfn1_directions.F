      submodule (wfn1) wfn1_directions
        integer :: ndir     ! the number of search directions
        integer :: ncurdir  ! the current number of directions
        integer :: nmo      ! the number of molecular orbitals
        integer :: nproc    ! the number processors
        integer :: iproc    ! the current processor rank
        logical :: oroot    ! is this the root processor?
        integer :: g_dirs   ! GA for search directions
        integer :: g_vec    ! GA for a single vector
        integer :: g_norms  ! the norms of the search directions
!
        double precision :: dnorm ! the vector norm
      contains
!
!     ------------------------------------------------------------------
!
      module subroutine wfn1_dir_init(wfn1_param)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      type(wfn1_prm), intent(in) :: wfn1_param
      double precision  :: done
      character(len=13) :: pname
      integer :: ilo, ihi, jlo, jhi, ii
      integer :: nea, neb
      pname   = "wfn1_dir_init"
      nmo     = wfn1_param%nmo
      nea     = wfn1_param%nea
      neb     = wfn1_param%neb
      if (nea.eq.neb) then
        ndir  = min(wfn1_param%ndir,max(0,2*nmo*(nmo-1)/2-24))
      else
        ndir  = min(wfn1_param%ndir,max(0,4*nmo*(nmo-1)/2-24))
      endif
      ncurdir = 0
      nproc   = ga_nnodes()
      iproc   = ga_nodeid()
      oroot   = iproc.eq.0
      done    = 1.0d0
!
!     Allocate the relevant global arrays
!
      if (.not.ga_create(MT_DBL,4*nmo**2,max(1,ndir),
     &                   "g_dirs",-1,ndir,g_dirs)) then
        call errquit(pname//" allocate g_dirs failed",
     &               wfn1_ibd_proc(4*(nmo**2)*ndir,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,4*nmo**2,1,
     &                   "g_vec",-1,1,g_vec)) then
        call errquit(pname//" allocate g_vec failed",
     &               wfn1_ibd_proc(4*(nmo**2),nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,1,max(1,ndir),
     &                   "g_norms",-1,-1,g_norms)) then
        call errquit(pname//" allocate g_norms failed",
     &               wfn1_ibd_proc(ndir,nproc),GA_ERR)
      endif
      call ga_zero(g_norms)
!
!     Set up unit vectors in g_dirs so we can invoke a simple
!     orthonormalization on it
!
      call ga_zero(g_dirs)
      call ga_distribution(g_dirs,iproc,ilo,ihi,jlo,jhi)
      do ii = max(ilo,jlo), min(ihi,jhi)
        call ga_put(g_dirs,ii,ii,ii,ii,done,1)
      enddo
      call ga_sync()
!
      end subroutine wfn1_dir_init
!
!     ------------------------------------------------------------------
!
      module subroutine wfn1_dir_orthog_vec(g_momo1,g_momo2,
     &                                      g_momo3,g_momo4)
      implicit none
#include "global.fh"
      integer, intent(in) :: g_momo1
      integer, intent(in) :: g_momo2
      integer, intent(in) :: g_momo3
      integer, intent(in) :: g_momo4
      integer :: nmo2
      nmo2 = nmo*nmo
      call ga_copy_patch('n',g_momo1,1,nmo, 1,nmo,
     &                       g_vec,  1,nmo2,1,   1)
      call ga_copy_patch('n',g_momo2,1,     nmo,   1,nmo,
     &                       g_vec,  nmo2+1,2*nmo2,1,  1)
      call ga_copy_patch('n',g_momo3,1,       nmo,   1,nmo,
     &                       g_vec,  2*nmo2+1,3*nmo2,1,  1)
      call ga_copy_patch('n',g_momo4,1,       nmo,   1,nmo,
     &                       g_vec,  3*nmo2+1,4*nmo2,1,  1)
      call wfn1_orthog_vec(g_dirs,ncurdir,g_vec,1)
      dnorm = ga_ddot(g_vec,g_vec)
      dnorm = sqrt(dnorm)
      call ga_copy_patch('n',g_vec,  1,nmo2,1,   1,
     &                       g_momo1,1,nmo, 1,nmo)
      call ga_copy_patch('n',g_vec,  nmo2+1,2*nmo2,1,  1,
     &                       g_momo2,1,     nmo,   1,nmo)
      call ga_copy_patch('n',g_vec,  2*nmo2+1,3*nmo2,1,  1,
     &                       g_momo3,1,       nmo,   1,nmo)
      call ga_copy_patch('n',g_vec,  3*nmo2+1,4*nmo2,1,  1,
     &                       g_momo4,1,       nmo,   1,nmo)
      end subroutine wfn1_dir_orthog_vec
!
!     ------------------------------------------------------------------
!
      module subroutine wfn1_dir_store_vec()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
      integer :: nmo4
      integer :: jj
      integer :: jjmax
      integer :: idummy
      double precision :: dnrmmax
      double precision, allocatable :: dnorms(:)
      idummy = 0
      nmo4 = 4*nmo*nmo
      if (ndir.le.0) return
      if (ncurdir.lt.ndir) then
!       Store the vector in addition to the already stored directions
        ncurdir = ncurdir + 1
        if (oroot) call ga_put(g_norms,1,1,ncurdir,ncurdir,dnorm,1)
        call ga_copy_patch('n',g_vec,1,nmo4,1,1,
     &                         g_dirs,1,nmo4,ncurdir,ncurdir)
      else
!       Find the vector with the largest norm to be replaced with the
!       new vector (larger norms correspond to likely worse harmonic
!       approximations).
        if (oroot) then
          allocate(dnorms(1:ndir))
          call ga_get(g_norms,1,1,1,ndir,dnorms,1)
          jjmax = 0
          dnrmmax = 0.0d0
          do jj = 1, ndir
            if (dnorms(jj).gt.dnrmmax) then
              jjmax = jj
              dnrmmax = dnorms(jj)
            endif
          enddo ! jj
          deallocate(dnorms)
          call ga_put(g_norms,1,1,jjmax,jjmax,dnorm,1)
        endif
        call ga_brdcst(001,jjmax,ma_sizeof(MT_BYTE,1,MT_INT),0)
        call ga_copy_patch('n',g_vec,1,nmo4,1,1,
     &                         g_dirs,1,nmo4,jjmax,jjmax)
      endif
      call wfn1_orthog_vecs(g_dirs,idummy,.false.)
      end subroutine wfn1_dir_store_vec
!
!     ------------------------------------------------------------------
!
      module subroutine wfn1_dir_final()
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      character(len=14) :: pname
      pname = "wfn1_dir_final"
      if (.not.ga_destroy(g_dirs)) then
        call errquit(pname//" deallocate g_dirs failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_vec)) then
        call errquit(pname//" deallocate g_vec failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(g_norms)) then
        call errquit(pname//" deallocate g_norms failed",0,GA_ERR)
      endif
      end subroutine wfn1_dir_final
!
!     ------------------------------------------------------------------
!
      end submodule wfn1_directions
