      submodule (wfn1) wfn1_gss_mvc
      implicit none
      contains
c-----------------------------------------------------------------------
c
!>
!> \brief Load the guess vectors from the movecs files
!>
!> In this routine we load the wave function from movecs files.
!> Because we have both natural orbitals and correlation functions
!> we expect there to be two files:
!> - one ending in ".movecs" which contains the natural orbitals
!> - one ending in ".covecs" which contains the correlation functions
!> If no ".covecs" file is available we will use the default guess
!> for the correlation functions.
!>
      module subroutine wfn1_guess_movecs(wfn1_param,wfn1_wave)
      implicit none
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
#include "errquit.fh"
#include "stdio.fh"
c
!>    The calculation parameters
      type(wfn1_prm), intent(in)    :: wfn1_param
!>    The wave function
      type(wfn1_wfn), intent(inout) :: wfn1_wave
c
      character*256 :: movecs_in
      character*256 :: covecs_in
      character*255 :: title_vecs
      character*255 :: bas_vecs
      character*4   :: file_scftype
      integer       :: nfilename
      integer       :: nbf_vecs
      integer       :: nsets
      integer       :: nmo_vecs(2)
      logical       :: status
c
      double precision, allocatable :: occ(:)
      double precision, allocatable :: eval(:)
c
      integer :: ii, nn ! counters
c
      logical, external :: movecs_read_header
      logical, external :: movecs_read
c
      allocate(occ(1:wfn1_param%nbf))
      allocate(eval(1:wfn1_param%nbf))
c
      nfilename = len(trim(wfn1_param%movecs_in))
      if (wfn1_param%movecs_in(nfilename-5:nfilename)=="movecs") then
        movecs_in = wfn1_param%movecs_in
        covecs_in = wfn1_param%movecs_in(1:nfilename-6)//"covecs"
      else
        write(luout,*)'movecs_in = ',trim(wfn1_param%movecs_in)
        call errquit("wfn1_guess_movecs: bad filenames",0,UERR)
      endif
c
c     == retrieve the natural orbitals ==
c
      status = movecs_read_header(movecs_in, title_vecs, bas_vecs,
     &         file_scftype, nbf_vecs, nsets, nmo_vecs, 2)
      if (status) then
        if (wfn1_param%nbf.eq.nbf_vecs.and.
     &      wfn1_param%nmo.eq.nmo_vecs(1)) then
          if (.not.movecs_read(movecs_in, 1, occ, eval,
     &        wfn1_wave%g_na)) then
            call errquit('wfn1_guess_movecs: movecs_read Na failed',
     &                   0,DISK_ERR)
          endif
          if (file_scftype.eq."UHF".or.
     &        file_scftype.eq."WFN1") then
            if (.not.movecs_read(movecs_in, 2, occ, eval,
     &          wfn1_wave%g_nb)) then
              call errquit('wfn1_guess_movecs: movecs_read Nb failed',
     &                     0,DISK_ERR)
            endif
          else
            call ga_copy(wfn1_wave%g_na,wfn1_wave%g_nb)
          endif
        endif
      endif
c
c     == retrieve or generate the correlation functions ==
c
      status = movecs_read_header(covecs_in, title_vecs, bas_vecs,
     &         file_scftype, nbf_vecs, nsets, nmo_vecs, 2)
      if (status) then
        if (wfn1_param%nbf.eq.nbf_vecs.and.
     &      wfn1_param%nmo.eq.nmo_vecs(1)) then
          if (.not.movecs_read(covecs_in, 1, occ, eval,
     &        wfn1_wave%g_ca)) then
            call errquit('wfn1_guess_movecs: movecs_read Ca failed',
     &                   0,DISK_ERR)
          else
            nn = 0
            do ii = 1, wfn1_param%nmo
              if (occ(ii).gt.0.0d0) then
                nn = nn + 1
                if (nn.gt.wfn1_param%nea) then
                  call errquit('wfn1_guess_movecs: inconsistent no. '//
     &                         'of electrons',nn,UERR)
                endif
                wfn1_wave%icnfa(nn) = ii
                wfn1_wave%dcnta(nn) = occ(ii)
              endif
            enddo
          endif
          if (.not.movecs_read(covecs_in, 2, occ, eval,
     &        wfn1_wave%g_cb)) then
            call errquit('wfn1_guess_movecs: movecs_read Cb failed',
     &                   0,DISK_ERR)
          else
            nn = 0
            do ii = 1, wfn1_param%nmo
              if (occ(ii).gt.0.0d0) then
                nn = nn + 1
                if (nn.gt.wfn1_param%neb) then
                  call errquit('wfn1_guess_movecs: inconsistent no. '//
     &                         'of electrons',nn,UERR)
                endif
                wfn1_wave%icnfb(nn) = ii
                wfn1_wave%dcntb(nn) = occ(ii)
              endif
            enddo
          endif
        endif
      else
        call wfn1_guess_corr(wfn1_param,wfn1_wave,eval,wfn1_param%nbf)
      endif
c
      deallocate(occ,eval)
      end
c
c-----------------------------------------------------------------------
c
      end submodule wfn1_gss_mvc
