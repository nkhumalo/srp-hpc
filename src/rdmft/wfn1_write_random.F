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
      complex(kind=8),  allocatable :: zvec(:)
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
      if (itp_gz1.ne.MT_DCPL) then
        call errquit("wfn1_write_random: gz_momo1 not complex",0,UERR)
      endif
      if (itp_gz2.ne.MT_DCPL) then
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
