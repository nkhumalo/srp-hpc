c-----------------------------------------------------------------------
c
      subroutine wfn1_print_analysis(wfn1_param,wfn1_inst,wfn1_wave,
     &                               en1a,en1b,ec1a,ec1b,g_2aa,g_2bb)
      implicit none
C>
C> \brief Print the wavefunction analysis
C>
#include "mafdecls.fh"
#include "stdio.fh"
#include "global.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
c
      type(wfn1_prm),     intent(in) :: wfn1_param
      type(wfn1_prminst), intent(in) :: wfn1_inst
      type(wfn1_wfn),     intent(in) :: wfn1_wave
!>    The \f$\alpha\f$ natural orbital energies
      double precision, intent(inout) :: en1a(1:wfn1_param%nmo)
!>    The \f$\alpha\f$ natural orbital energies
      double precision, intent(inout) :: en1b(1:wfn1_param%nmo)
!>    The \f$\alpha\f$ correlation function energies
      double precision, intent(inout) :: ec1a(1:wfn1_param%nmo)
!>    The \f$\alpha\f$ correlation function energies
      double precision, intent(inout) :: ec1b(1:wfn1_param%nmo)
!>    The \f$\alpha\alpha\f$ 2-matrix occupation numbers
      integer,          intent(inout) :: g_2aa
!>    The \f$\beta\beta\f$ 2-matrix occupation numbers
      integer,          intent(inout) :: g_2bb
!>    The natural orbital occupation numbers
      double precision, allocatable   :: ocna(:), ocnb(:)
      double precision, allocatable   :: occa(:), occb(:)
      double precision, allocatable   :: ci(:), cj(:)
c
      double precision :: tanalyze
      double precision :: valmin
      integer          :: imin
      integer          :: iproc ! current rank
      integer          :: nproc ! number of all ranks
c
      integer, allocatable :: irrep(:)
c
      integer :: ii,  jj  ! counters
      integer :: ilo, ihi ! range limits
      integer :: nmo      ! local copy of number of orbitals
      integer :: ne_a     ! local copy of number of alpha electrons
      integer :: ne_b     ! local copy of number of alpha electrons
c
      character*255 blob
c
      iproc    = ga_nodeid()
      nproc    = ga_nnodes()
      tanalyze = 0.025d0
      ilo  = 1
      ihi  = wfn1_param%nmo
      nmo  = wfn1_param%nmo
      ne_a = wfn1_param%nea
      ne_b = wfn1_param%neb
c
      allocate(ocna(1:nmo),ocnb(1:nmo))
      allocate(occa(1:nmo),occb(1:nmo))
      allocate(ci(1:nmo),cj(1:nmo))
      allocate(irrep(1:nmo))
      irrep = 0
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,wfn1_param%nea,
     &                   wfn1_wave%icnfa,wfn1_wave%dcnta,ocna)
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,wfn1_param%neb,
     &                   wfn1_wave%icnfb,wfn1_wave%dcntb,ocnb)
c
c     2-matrix occupation numbers
c
      call wfn1_calc_ocx(wfn1_wave%g_ca,nmo,ne_a,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,g_2aa)
      call wfn1_calc_ocx(wfn1_wave%g_cb,nmo,ne_b,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,g_2bb)
      call wfn1_add_occd(nmo,ocna,g_2aa)
      call wfn1_add_occd(nmo,ocnb,g_2bb)
      call ga_sync()
c
c     call wfn1_print_occ(nmo,occa,occb)
c
      if (iproc.eq.0) then
        blob = "WFN1 Final Natural orbital energies and occupations"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        write(LuOut,'(12x,"  Alpha ",4x,"  Beta  ",8x,
     &                    "  Alpha ",4x,"  Beta  ")')
        write(LuOut,'(12x,"--------",4x,"--------",8x,
     &                    "--------",4x,"--------")')
        do ii = 1, nmo
          write(LuOut,'(i8,f12.4,f12.4,4x,f12.4,f12.4)')ii,
     &          en1a(ii),en1b(ii),ocna(ii),ocnb(ii)
        enddo
cDEBUG
        do ii = 1, nmo
          write(80,'(2x,f14.6,$)')ocna(ii)
        enddo
        do ii = 1, nmo
          write(80,'(2x,f14.6,$)')ocnb(ii)
        enddo
        close(80)
cDEBUG
      endif
c
c     Sort correlation functions by energy
c
c     Alpha
c
      occa = 0.0d0
      occb = 0.0d0
      do ii = 1, wfn1_param%nea
        occa(wfn1_wave%icnfa(ii)) = wfn1_wave%dcnta(ii)
      enddo
      do ii = 1, wfn1_param%neb
        occb(wfn1_wave%icnfb(ii)) = wfn1_wave%dcntb(ii)
      enddo
      if (iproc.eq.0) then
        do ii = 1, nmo-1
          valmin = ec1a(ii)
          imin   = ii
          do jj = ii+1, nmo
            if (ec1a(jj).lt.valmin) then
              valmin = ec1a(jj)
              imin   = jj
            endif
          enddo
          if (imin.ne.ii) then
            call ga_get(wfn1_wave%g_ca,1,nmo,ii,ii,ci,nmo)
            call ga_get(wfn1_wave%g_ca,1,nmo,imin,imin,cj,nmo)
            call ga_put(wfn1_wave%g_ca,1,nmo,ii,ii,cj,nmo)
            call ga_put(wfn1_wave%g_ca,1,nmo,imin,imin,ci,nmo)
            ec1a(imin) = ec1a(ii)
            ec1a(ii)   = valmin
            valmin     = occa(imin)
            occa(imin) = occa(ii)
            occa(ii)   = valmin
          endif
        enddo
      endif
c
c     Beta
c
      if (iproc.eq.nproc-1) then
        do ii = 1, nmo-1
          valmin = ec1b(ii)
          imin   = ii
          do jj = ii+1, nmo
            if (ec1b(jj).lt.valmin) then
              valmin = ec1b(jj)
              imin   = jj
            endif
          enddo
          if (imin.ne.ii) then
            call ga_get(wfn1_wave%g_cb,1,nmo,ii,ii,ci,nmo)
            call ga_get(wfn1_wave%g_cb,1,nmo,imin,imin,cj,nmo)
            call ga_put(wfn1_wave%g_cb,1,nmo,ii,ii,cj,nmo)
            call ga_put(wfn1_wave%g_cb,1,nmo,imin,imin,ci,nmo)
            ec1b(imin) = ec1b(ii)
            ec1b(ii)   = valmin
            valmin     = occb(imin)
            occb(imin) = occb(ii)
            occb(ii)   = valmin
          endif
        enddo
      endif
      call ga_sync()
      call ga_brdcst(1100,ec1a,ma_sizeof(MT_DBL,nmo,MT_BYTE),0)
      call ga_brdcst(1101,occa,ma_sizeof(MT_DBL,nmo,MT_BYTE),0)
      call ga_brdcst(1102,ec1b,ma_sizeof(MT_DBL,nmo,MT_BYTE),nproc-1)
      call ga_brdcst(1103,occb,ma_sizeof(MT_DBL,nmo,MT_BYTE),nproc-1)
c
      if (iproc.eq.0) then
        blob="WFN1 Final Correlation function energies and occupations"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        write(LuOut,'(12x,"  Alpha ",4x,"  Beta  ",8x,
     &                    "  Alpha ",4x,"  Beta  ")')
        write(LuOut,'(12x,"--------",4x,"--------",8x,
     &                    "--------",4x,"--------")')
        do ii = ilo, ihi
          if (abs(ec1a(ii)).gt.1.0d5.or.abs(ec1b(ii)).gt.1.0d5) then
            write(LuOut,'(i8,e12.1,e12.1,4x,f12.4,f12.4)')ii,
     &            ec1a(ii),ec1b(ii),occa(ii),occb(ii)
          else
            write(LuOut,'(i8,f12.4,f12.4,4x,f12.4,f12.4)')ii,
     &            ec1a(ii),ec1b(ii),occa(ii),occb(ii)
          endif
        enddo
      endif
c
      if (iproc.eq.0) then
        blob = "WFN1 Final alpha-alpha 2-matrix occupation numbers"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(g_2aa)
      if (iproc.eq.0) then
        call util_flush(LuOut)
      endif
c
      if (iproc.eq.0) then
        blob = "WFN1 Final beta-beta 2-matrix occupation numbers"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(g_2bb)
      if (iproc.eq.0) then
        call util_flush(LuOut)
      endif
c
      if (iproc.eq.0) then
        blob = "WFN1 Final alpha Correlation Function Analysis"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_ca)
      if (iproc.eq.0) then
        call util_flush(LuOut)
      endif
c
      if (iproc.eq.0) then
        blob = "WFN1 Final beta Correlation Function Analysis"
        write(LuOut,*)
        call util_print_centered(LuOut,blob, 40, .true.)
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_cb)
      if (iproc.eq.0) then
        call util_flush(LuOut)
      endif
c
      blob = "WFN1 Final Alpha Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_na,blob,.true.,en1a,.false.,
     &                       irrep,.true.,ocna)
c
      blob = "WFN1 Final Beta Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_nb,blob,.true.,en1b,.false.,
     &                       irrep,.true.,ocnb)
c
      deallocate(occa,occb)
      deallocate(ocna,ocnb)
      deallocate(ci,cj)
c
      end
c
c-----------------------------------------------------------------------
