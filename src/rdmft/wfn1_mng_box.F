!-----------------------------------------------------------------------
!>
!> \brief Allocate box global arrays
!>
      subroutine wfn1_allocate_box(wfn1_param,wfn1_box)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "mafdecls.fh"
!>    The wavefunction parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The Monte Carlo box
      type(wfn1_bx),  intent(inout) :: wfn1_box
!>
      integer :: nmo
      integer :: nproc
      character(len=18), parameter :: pname = "wfn1_allocate_box:"
!>
      nmo = wfn1_param%nmo
      nproc = ga_nnodes()
      if (.not.ga_create(MT_DBL,nmo,nmo,"box%g_na",-1,-1,wfn1_box%g_na))
     &then
        call errquit(pname//" allocate box%g_na failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,nmo,nmo,"box%g_nb",-1,-1,wfn1_box%g_nb))
     &then
        call errquit(pname//" allocate box%g_nb failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,nmo,nmo,"box%g_ca",-1,-1,wfn1_box%g_ca))
     &then
        call errquit(pname//" allocate box%g_ca failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,nmo,nmo,"box%g_cb",-1,-1,wfn1_box%g_cb))
     &then
        call errquit(pname//" allocate box%g_cb failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      end subroutine wfn1_allocate_box
!
!-----------------------------------------------------------------------
!>
!> \brief Deallocate box global arrays
!>
      subroutine wfn1_deallocate_box(wfn1_box)
      implicit none
#include "global.fh"
#include "errquit.fh"
!>    The Monte Carlo box
      type(wfn1_bx),  intent(inout) :: wfn1_box
!>
      character(len=20), parameter :: pname = "wfn1_deallocate_box:"
!>
      if (.not.ga_destroy(wfn1_box%g_na)) then
        call errquit(pname//" deallocate box%g_na failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_box%g_nb)) then
        call errquit(pname//" deallocate box%g_nb failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_box%g_ca)) then
        call errquit(pname//" deallocate box%g_ca failed",0,GA_ERR)
      endif
      if (.not.ga_destroy(wfn1_box%g_cb)) then
        call errquit(pname//" deallocate box%g_cb failed",0,GA_ERR)
      endif
      end subroutine wfn1_deallocate_box
!
!-----------------------------------------------------------------------
!>
!> \brief Box initialize to a specified value
!>
      subroutine wfn1_init_box(box,val)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
      target :: dbl_mb
!>    The Monte Carlo box
      type(wfn1_bx), intent(inout) :: box
!>    The value to initialize to
      real(kind=8),  intent(in)    :: val
!
      integer :: ilo, ihi, jlo, jhi
      integer :: nsize, ldi, ldj, idx
      integer :: iproc
      real(kind=8), pointer :: ptr(:,:)
!
      iproc = ga_nodeid()
!
      call ga_distribution(box%g_na,iproc,ilo,ihi,jlo,jhi)
      ldi = ihi-ilo+1
      ldj = jhi-jlo+1
      nsize = ldi*ldj
      call ga_access(box%g_na,ilo,ihi,jlo,jhi,idx,ldi)
      ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsize-1)
      ptr = val
!
      call ga_distribution(box%g_nb,iproc,ilo,ihi,jlo,jhi)
      ldi = ihi-ilo+1
      ldj = jhi-jlo+1
      nsize = ldi*ldj
      call ga_access(box%g_nb,ilo,ihi,jlo,jhi,idx,ldi)
      ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsize-1)
      ptr = val
!
      call ga_distribution(box%g_ca,iproc,ilo,ihi,jlo,jhi)
      ldi = ihi-ilo+1
      ldj = jhi-jlo+1
      nsize = ldi*ldj
      call ga_access(box%g_ca,ilo,ihi,jlo,jhi,idx,ldi)
      ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsize-1)
      ptr = val
!
      call ga_distribution(box%g_cb,iproc,ilo,ihi,jlo,jhi)
      ldi = ihi-ilo+1
      ldj = jhi-jlo+1
      nsize = ldi*ldj
      call ga_access(box%g_cb,ilo,ihi,jlo,jhi,idx,ldi)
      ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsize-1)
      ptr = val
!
      end subroutine wfn1_init_box
!
!-----------------------------------------------------------------------
!
!> \brief Find the maximum of every dimension in box
!>
      subroutine wfn1_max_box(wfn1_box,g_rna,g_rnb,g_rca,g_rcb)
      implicit none
!>    The Monte Carlo box
      type(wfn1_bx), intent(inout) :: wfn1_box
!>    The alpha natural orbital rotation matrix
      integer,       intent(in)    :: g_rna
!>    The beta  natural orbital rotation matrix
      integer,       intent(in)    :: g_rnb
!>    The alpha correlation function rotation matrix
      integer,       intent(in)    :: g_rca
!>    The beta  correlation function rotation matrix
      integer,       intent(in)    :: g_rcb
      call wfn1_max_arrays(wfn1_box%g_na,g_rna)
      call wfn1_max_arrays(wfn1_box%g_nb,g_rnb)
      call wfn1_max_arrays(wfn1_box%g_ca,g_rca)
      call wfn1_max_arrays(wfn1_box%g_cb,g_rcb)
      end subroutine wfn1_max_box
!
!-----------------------------------------------------------------------
!
!> \brief Find the maximum of arrays
!>
      subroutine wfn1_max_arrays(g_max,g_a)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      target dbl_mb
!>    The array that keeps the maximum values
      integer, intent(inout) :: g_max
!>    The arrray with new values
      integer, intent(in)    :: g_a
!
      character(len=16), parameter :: pname = "wfn1_max_arrays"
      integer itype, ilo, ihi, jlo, jhi, iln, jln
      integer idx,   ldx, nsz, ii,  jj
      integer iproc
      real(kind=8), pointer :: ptr(:,:)
!
      real(kind=8),    allocatable :: buf(:,:)
      complex(kind=8), allocatable :: zbuf(:,:)
!
      iproc = ga_nodeid()
      call ga_distribution(g_max,iproc,ilo,ihi,jlo,jhi)
      call ga_inquire(g_a,itype,iln,jln)
      if      (itype.eq.MT_DBL) then
        allocate(buf(ilo:ihi,jlo:jhi))
        call ga_get(g_a,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
        call ga_access(g_max,ilo,ihi,jlo,jhi,idx,ldx)
        nsz = (ihi-ilo+1)*(jhi-jlo+1)
        ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsz-1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            ptr(ii,jj) = max(ptr(ii,jj),abs(buf(ii,jj)))
          enddo
        enddo
        deallocate(buf)
      else if (itype.eq.MT_DCPL) then
        allocate(zbuf(ilo:ihi,jlo:jhi))
        call ga_get(g_a,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
        call ga_access(g_max,ilo,ihi,jlo,jhi,idx,ldx)
        nsz = (ihi-ilo+1)*(jhi-jlo+1)
        ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsz-1)
        do jj = jlo, jhi
          do ii = ilo, ihi
            ptr(ii,jj) = max(ptr(ii,jj),abs(zbuf(ii,jj)))
          enddo
        enddo
        deallocate(zbuf)
      else
        call errquit(pname//" unexpected type",itype,UERR)
      endif
      end subroutine wfn1_max_arrays
!
!-----------------------------------------------------------------------
!
!> \brief calculate the volume of an array of the Monte Carlo box
!>
      function wfn1_vol_array(g_a) result(volume)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      target :: dbl_mb
!>    The Monte Carlo box
      integer, intent(in) :: g_a
!>    The box volume
      real(kind=8) :: volume
!
      integer :: nproc
      integer :: iproc
      integer :: ilo, ihi, jlo, jhi, idx, ldx
      integer :: nsz, ii, jj
!
      real(kind=8), allocatable :: buf(:)
      real(kind=8), pointer     :: ptr(:,:)
!
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      allocate(buf(0:nproc-1))
!
      call ga_distribution(g_a,iproc,ilo,ihi,jlo,jhi)
      call ga_access(g_a,ilo,ihi,jlo,jhi,idx,ldx)
      nsz = (ihi-ilo+1)*(jhi-jlo+1)
      ptr(ilo:ihi,jlo:jhi) => dbl_mb(idx:idx+nsz-1)
      volume = 1.0d0
      do jj = jlo, jhi
        do ii = max(ilo,jj+1), ihi
          volume = volume * ptr(ii,jj)
        enddo
      enddo
      buf = 0.0d0
      buf(iproc) = volume
      call ga_dgop(2000,buf,nproc,"+")
      volume = 1.0d0
      do ii = 0, nproc-1
        volume = volume * buf(ii)
      enddo
!
      deallocate(buf)
      end function wfn1_vol_array
!
!-----------------------------------------------------------------------
!
!> \brief calculate the volume of the Monte Carlo box
!>
      function wfn1_vol_box(wfn1_box) result(volume)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The Monte Carlo box
      type(wfn1_bx), intent(in) :: wfn1_box
!>    The box volume
      real(kind=8) :: volume
!
      volume = wfn1_vol_array(wfn1_box%g_na)
      volume = volume * wfn1_vol_array(wfn1_box%g_nb)
      volume = volume * wfn1_vol_array(wfn1_box%g_ca)
      volume = volume * wfn1_vol_array(wfn1_box%g_cb)
      end function wfn1_vol_box
!
!-----------------------------------------------------------------------
!
!> \brief normalize box volume
!>
      subroutine wfn1_normalize_box(wfn1_param,wfn1_box)
      implicit none
!>    Wavefunction parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The Monte Carlo box
      type(wfn1_bx),  intent(inout) :: wfn1_box
!
      integer :: numdims
      integer :: nmo
      real(kind=8) :: volume
      real(kind=8) :: scalef
!
      nmo = wfn1_param%nmo
      numdims = 4*(nmo*(nmo-1)/2)
      volume = wfn1_vol_box(wfn1_box)
      scalef = (1.0d0/volume)**(1.0d0/dble(numdims))
      call ga_scale(wfn1_box%g_na,scalef)
      call ga_scale(wfn1_box%g_nb,scalef)
      call ga_scale(wfn1_box%g_ca,scalef)
      call ga_scale(wfn1_box%g_cb,scalef)
      end subroutine wfn1_normalize_box
!
!-----------------------------------------------------------------------
!
!> \brief Change the Monte Carlo box shape
!>
!> To change the box shape we look at the length of the sides of the 
!> rejected steps that exceed the length of the accepted step. From this
!> we assume that the size of the box in those dimensions are too big
!> and shrink them. At the end we renormalize the box back to unit
!> volume.
!>
      subroutine wfn1_reshape_box(wfn1_param,wfn1_box,wfn1_rejects,
     &                            wfn1_accept)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The wavefunction parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The Monte Carlo box
      type(wfn1_bx),  intent(inout) :: wfn1_box
!>    The box outline of all rejected steps
      type(wfn1_bx),  intent(in)    :: wfn1_rejects
!>    The box outline of the accepted step
      type(wfn1_bx),  intent(in)    :: wfn1_accept
!
      call wfn1_reshape_array(wfn1_box%g_na,wfn1_rejects%g_na,
     &                        wfn1_accept%g_na)
      call wfn1_reshape_array(wfn1_box%g_nb,wfn1_rejects%g_nb,
     &                        wfn1_accept%g_nb)
      call wfn1_reshape_array(wfn1_box%g_ca,wfn1_rejects%g_ca,
     &                        wfn1_accept%g_ca)
      call wfn1_reshape_array(wfn1_box%g_cb,wfn1_rejects%g_cb,
     &                        wfn1_accept%g_cb)
      call wfn1_normalize_box(wfn1_param,wfn1_box)
!
      end subroutine wfn1_reshape_box
!
!-----------------------------------------------------------------------
!
      subroutine wfn1_reshape_array(g_box,g_rejects,g_accept)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      target :: dbl_mb
!>    The Monte Carlo array
      integer, intent(inout) :: g_box
!>    The array outline of all rejected steps
      integer, intent(in)    :: g_rejects
!>    The array outline of the accepted step
      integer, intent(in)    :: g_accept
!
      integer :: iproc
      integer :: iblo, ibhi, jblo, jbhi
      integer :: irlo, irhi, jrlo, jrhi
      integer :: ialo, iahi, jalo, jahi
      integer :: nsz, idb, ldb, idr, ldr, ida, lda
      integer :: ii, jj
      real(kind=8), pointer :: pb(:,:)
      real(kind=8), pointer :: pr(:,:)
      real(kind=8), pointer :: pa(:,:)
!
      iproc = ga_nodeid()
      call ga_distribution(g_box,iproc,iblo,ibhi,jblo,jbhi)
      call ga_distribution(g_rejects,iproc,irlo,irhi,jrlo,jrhi)
      call ga_distribution(g_accept,iproc,ialo,iahi,jalo,jahi)
      if (iblo.ne.irlo.or.iblo.ne.ialo.or.
     &    ibhi.ne.irhi.or.ibhi.ne.iahi.or.
     &    jblo.ne.jrlo.or.jblo.ne.jalo.or.
     &    jbhi.ne.jrhi.or.jbhi.ne.jahi) then
        call errquit("distributions mismatch",0,UERR)
      endif
      nsz = (ibhi-iblo+1)*(jbhi-jblo+1)
      call ga_access(g_box,iblo,ibhi,jblo,jbhi,idb,ldb)
      pb(iblo:ibhi,jblo:jbhi) => dbl_mb(idb:idb+nsz-1)
      call ga_access(g_rejects,irlo,irhi,jrlo,jrhi,idr,ldr)
      pr(irlo:irhi,jrlo:jrhi) => dbl_mb(idr:idr+nsz-1)
      call ga_access(g_accept,ialo,iahi,jalo,jahi,ida,lda)
      pa(ialo:iahi,jalo:jahi) => dbl_mb(ida:ida+nsz-1)
      do jj = jblo, jbhi
        do ii = iblo, ibhi
          if (pr(ii,jj).gt.pa(ii,jj)) then
            pb(ii,jj) = pb(ii,jj)*sqrt(pa(ii,jj)/pr(ii,jj))
          endif
        enddo
      enddo
      end subroutine wfn1_reshape_array
!
!-----------------------------------------------------------------------
!
!> \brief Change the Monte Carlo box shape towards square
!>
!> Because the energy is invariant to some rotations the box can be
!> skewed to be exceedingly large in some dimensions and small in
!> others. Resquaring the box might help.
!>
      subroutine wfn1_resquare_box(wfn1_param,wfn1_box)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The wavefunction parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The Monte Carlo box
      type(wfn1_bx),  intent(inout) :: wfn1_box
!
      call wfn1_resquare_array(wfn1_box%g_na)
      call wfn1_resquare_array(wfn1_box%g_nb)
      call wfn1_resquare_array(wfn1_box%g_ca)
      call wfn1_resquare_array(wfn1_box%g_cb)
      call wfn1_normalize_box(wfn1_param,wfn1_box)
!
      end subroutine wfn1_resquare_box
!
!-----------------------------------------------------------------------
!
      subroutine wfn1_resquare_array(g_box)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      target :: dbl_mb
!>    The Monte Carlo array
      integer, intent(inout) :: g_box
!
      integer :: iproc
      integer :: iblo, ibhi, jblo, jbhi
      integer :: nsz, idb, ldb, idr, ldr, ida, lda
      integer :: ii, jj
      real(kind=8), pointer :: pb(:,:)
!
      iproc = ga_nodeid()
      call ga_distribution(g_box,iproc,iblo,ibhi,jblo,jbhi)
      nsz = (ibhi-iblo+1)*(jbhi-jblo+1)
      call ga_access(g_box,iblo,ibhi,jblo,jbhi,idb,ldb)
      pb(iblo:ibhi,jblo:jbhi) => dbl_mb(idb:idb+nsz-1)
      do jj = jblo, jbhi
        do ii = iblo, ibhi
          pb(ii,jj) = sqrt(pb(ii,jj))
        enddo
      enddo
      end subroutine wfn1_resquare_array
!
!-----------------------------------------------------------------------
