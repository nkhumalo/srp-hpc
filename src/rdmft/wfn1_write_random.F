!-----------------------------------------------------------------------
!>
!> \brief Write random box to file
!>
!> Write the random box to a file so we can use the data somewhere else.
!>
      subroutine wfn1_write_random(filename,no,g_momo1,g_momo2,
     &                                         gz_momo1,gz_momo2)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The filename prefix
      character(len=*), intent(in) :: filename
!>    The Fock build number
      integer,          intent(in) :: no
!>    The alpha natural orbital random box
      integer,          intent(in) :: g_momo1
!>    The beta  natural orbital random box
      integer,          intent(in) :: g_momo2
!>    The alpha correlation function random box
      integer,          intent(in) :: gz_momo1
!>    The beta  correlation function random box
      integer,          intent(in) :: gz_momo2
!
!     Local
!
      integer :: itp_g1, itp_g2, itp_gz1, itp_gz2
      integer :: idm_g1, idm_g2, idm_gz1, idm_gz2
      integer :: jdm_g1, jdm_g2, jdm_gz1, jdm_gz2
!
      integer :: ii, jj, nmo
      integer :: iunit
!
      real(kind=8),     allocatable :: dvec(:)
c     complex(kind=8),  allocatable :: zvec(:)
      real   (kind=8),  allocatable :: zvec(:)
      character(len=:), allocatable :: fname
!
      iunit = 10
      call ga_inquire(g_momo1, itp_g1, idm_g1, jdm_g1)
      call ga_inquire(g_momo2, itp_g2, idm_g2, jdm_g2)
      call ga_inquire(gz_momo1,itp_gz1,idm_gz1,jdm_gz1)
      call ga_inquire(gz_momo2,itp_gz2,idm_gz2,jdm_gz2)
!
      if (itp_g1.ne.MT_DBL) then
        call errquit("wfn1_write_random: g_momo1 not double",0,UERR)
      endif
      if (itp_g2.ne.MT_DBL) then
        call errquit("wfn1_write_random: g_momo2 not double",0,UERR)
      endif
      if (itp_gz1.ne.MT_DBL) then
        call errquit("wfn1_write_random: gz_momo1 not complex",0,UERR)
      endif
      if (itp_gz2.ne.MT_DBL) then
        call errquit("wfn1_write_random: gz_momo2 not complex",0,UERR)
      endif
      if (idm_g1.eq.idm_g2.and.idm_g1.eq.idm_gz1.and.idm_g1.eq.idm_gz2)
     &then
        nmo = idm_g1
      else
        call errquit("wfn1_write_random: inconsistent dimensions",
     &               idm_g1,UERR)
      endif
      if (idm_g1 .ne.jdm_g1 .or.idm_g2 .ne.jdm_g2 .or.
     &    idm_gz1.ne.jdm_gz1.or.idm_gz2.ne.jdm_gz2) then
        call errquit("wfn1_write_random: non-square matrices",0,UERR)
      endif
!
      call ga_sync()
      if (ga_nodeid().eq.0) then
        ii = len(trim(filename))+6
        allocate(character(len=ii) :: fname)
        write(fname,'(a,".",i0.5)')filename,no
        open(unit=iunit,file=fname)
        allocate(dvec(1:nmo),zvec(1:nmo))
        write(iunit,*)nmo
!
        do jj = 1, nmo
          call ga_get(g_momo1,1,nmo,jj,jj,dvec,nmo)
          do ii = 1,nmo
            write(iunit,*)dvec(ii)
          enddo
        enddo
        do jj = 1, nmo
          call ga_get(g_momo2,1,nmo,jj,jj,dvec,nmo)
          do ii = 1,nmo
            write(iunit,*)dvec(ii)
          enddo
        enddo
        do jj = 1, nmo
          call ga_get(gz_momo1,1,nmo,jj,jj,zvec,nmo)
          do ii = 1,nmo
            write(iunit,*)zvec(ii)
          enddo
        enddo
        do jj = 1, nmo
          call ga_get(gz_momo2,1,nmo,jj,jj,zvec,nmo)
          do ii = 1,nmo
            write(iunit,*)zvec(ii)
          enddo
        enddo
!
        close(unit=iunit)
        deallocate(dvec,zvec)
        deallocate(fname)
        call util_flush(6)
      endif
      call ga_sync()
      end subroutine wfn1_write_random
!>
!-----------------------------------------------------------------------
!>
!> \brief Write the parabolic approximation to a file in GNU Plot format
!>
!> Write the parabolic approximation to a file so we can visualize it
!> with GNU plot.
!>
      subroutine wfn1_write_parabola(filename,no,aa,bb,cc)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The filename prefix
      character(len=*), intent(in) :: filename
!>    The Fock build number
      integer,          intent(in) :: no
!>    The alpha natural orbital random box
      real(kind=8),     intent(in) :: aa
!>    The beta  natural orbital random box
      real(kind=8),     intent(in) :: bb
!>    The alpha correlation function random box
      real(kind=8),     intent(in) :: cc
!
!     Local
!
      integer :: iunit
      integer :: ii
!
      character(len=:), allocatable :: fname
!
      iunit = 10
!
      call ga_sync()
      if (ga_nodeid().eq.0) then
        ii = len(trim(filename))+6
        allocate(character(len=ii) :: fname)
        write(fname,'(a,".",i0.5)')filename,no
        open(unit=iunit,file=fname)
        write(iunit,*)"aa=",aa
        write(iunit,*)"bb=",bb
        write(iunit,*)"cc=",cc
        write(iunit,*)"p(x)=aa*x*x+bb*x+cc"
        close(unit=iunit)
        deallocate(fname)
        call util_flush(6)
      endif
      call ga_sync()
      end subroutine wfn1_write_parabola
!
!-----------------------------------------------------------------------
!>
!> \brief Write the known data points
!>
!> Write the known data points so we can visualize them with GNU Plot.
!>
      subroutine wfn1_write_points(filename,no,xm,x0,xp,xx,Em1,E0,Ep1,
     &                                                     Emin,Ewf)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The filename prefix
      character(len=*), intent(in) :: filename
!>    The Fock build number
      integer,          intent(in) :: no
!>    The minimum x
      real(kind=8),     intent(in) :: xm
!>    The zero x
      real(kind=8),     intent(in) :: x0
!>    The maximum x
      real(kind=8),     intent(in) :: xp
!>    The maximum x
      real(kind=8),     intent(in) :: xx
!>    The energy at xm
      real(kind=8),     intent(in) :: Em1
!>    The energy at x0
      real(kind=8),     intent(in) :: E0
!>    The energy at xp
      real(kind=8),     intent(in) :: Ep1
!>    The energy at xx
      real(kind=8),     intent(in) :: Emin, Ewf
!
!     Local
!
      integer :: iunit
      integer :: ii
!
      character(len=:), allocatable :: fname
!
      iunit = 10
!
      call ga_sync()
      if (ga_nodeid().eq.0) then
        ii = len(trim(filename))+6
        allocate(character(len=ii) :: fname)
        write(fname,'(a,".",i0.5)')filename,no
        open(unit=iunit,file=fname)
        write(iunit,*)xm," ",Em1
        if (xx.lt.x0) then
          write(iunit,*)xx," ",Emin
          write(iunit,*)xx," ",Ewf
        endif
        write(iunit,*)x0," ",E0
        if (x0.lt.xx) then
          write(iunit,*)xx," ",Emin
          write(iunit,*)xx," ",Ewf
        endif
        write(iunit,*)xp," ",Ep1
        close(unit=iunit)
        deallocate(fname)
        call util_flush(6)
      endif
      call ga_sync()
      end subroutine wfn1_write_points
!
!-----------------------------------------------------------------------
!>
!> \brief Write the data points along the line
!>
      subroutine wfn1_write_line(filename,no,
     &           wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &           wfn1_trial_wv,wfn1_trial_en,
     &           g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &           g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,g_momo5,g_momo6,
     &           gz_momo1,gz_momo2,gz_momo3,gz_momo4,
     &           g_4t1,g_4t2)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!>    The filename prefix
      character(len=*), intent(in) :: filename
!>    The Fock build number
      integer,          intent(in) :: no
!
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
      type(wfn1_wfn),     intent(inout) :: wfn1_trial_wv
      type(wfn1_en),      intent(out)   :: wfn1_trial_en
!
      integer, intent(in) :: g_bfbf1 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf2 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf3 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf4 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfmo  ! nbf x nmo scratch GA
      integer, intent(in) :: g_momo1 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo2 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo3 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo4 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo5 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo6 ! nmo x nmo scratch GA
      integer, intent(in) :: gz_momo1 ! nmo x nmo scratch GA
      integer, intent(in) :: gz_momo2 ! nmo x nmo scratch GA
      integer, intent(in) :: gz_momo3 ! nmo x nmo scratch GA
      integer, intent(in) :: gz_momo4 ! nmo x nmo scratch GA
      integer, intent(in) :: g_4t1   ! nbf x nbf x nbf x nbf scratch GA
      integer, intent(in) :: g_4t2   ! nbf x nbf x nbf x nbf scratch GA
!
!     Local
!
      integer             :: gz_momo5 ! nmo x nmo scratch GA
      integer             :: gz_momo6 ! nmo x nmo scratch GA
!
      integer :: iunit
      integer :: nmo
      real(kind=8) :: xx
      real(kind=8) :: Ewf_trial
      integer :: iinum,ii
      integer :: nproc
!
      character(len=15), parameter :: pname = "wfn1_write_line"
      character(len=:), allocatable :: fname
!
      iunit = 10
      nmo = wfn1_param%nmo
      iinum = 25
      nproc = ga_nnodes()
!
      call ga_sync()
      if (.not.ga_create(MT_DBL,nmo,nmo,"gz_momo5",-1,-1,gz_momo5)) then
        call errquit(pname//" allocate gz_momo5 failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_create(MT_DBL,nmo,nmo,"gz_momo6",-1,-1,gz_momo6)) then
        call errquit(pname//" allocate gz_momo6 failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (ga_nodeid().eq.0) then
        ii = len(trim(filename))+6
        allocate(character(len=ii) :: fname)
        write(fname,'(a,".",i0.5)')filename,no
        open(unit=iunit,file=fname)
      endif
      do ii = -iinum, iinum
        xx = dble(ii)/dble(iinum)
        call ga_copy(g_momo1,g_momo3)
        call ga_copy(g_momo2,g_momo4)
        call ga_scale(g_momo3,xx)
        call ga_scale(g_momo4,xx)
        call wfn1_skew_to_rot(g_momo3,g_momo5)
        call wfn1_skew_to_rot(g_momo4,g_momo6)
        call ga_copy(gz_momo1,gz_momo3)
        call ga_copy(gz_momo2,gz_momo4)
        call ga_scale(gz_momo3,xx)
        call ga_scale(gz_momo4,xx)
        call wfn1_skew_to_rot(gz_momo3,gz_momo5)
        call wfn1_skew_to_rot(gz_momo4,gz_momo6)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,
     &                 g_momo5,0.0d0,wfn1_trial_wv%g_na)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,
     &                 g_momo6,0.0d0,wfn1_trial_wv%g_nb)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &                 gz_momo5,0.0d0,wfn1_trial_wv%g_ca)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &                 gz_momo6,0.0d0,wfn1_trial_wv%g_cb)
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                wfn1_trial_en,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,g_momo3,g_momo4,g_momo5,g_momo6,
     &                g_4t1,g_4t2)
        Ewf_trial  = wfn1_trial_en%E_tot
        if (ga_nodeid().eq.0) then
          write(*,*)'HVD:',ii,"  ",xx,"  ",Ewf_trial
          call util_flush(6)
          write(iunit,*)xx,"  ",Ewf_trial
        endif
      enddo ! ii
      if (ga_nodeid().eq.0) then
        close(unit=iunit)
        deallocate(fname)
        call util_flush(6)
      endif
      if (.not.ga_destroy(gz_momo5)) then
        call errquit(pname//" deallocate gz_momo5 failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      if (.not.ga_destroy(gz_momo6)) then
        call errquit(pname//" deallocate gz_momo6 failed",
     &               wfn1_ibd_proc(nmo**2,nproc),GA_ERR)
      endif
      call ga_sync()
      end subroutine wfn1_write_line
!>
!-----------------------------------------------------------------------
