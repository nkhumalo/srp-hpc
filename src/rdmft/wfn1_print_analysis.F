c-----------------------------------------------------------------------
c
      subroutine wfn1_print_analysis(wfn1_param,wfn1_inst,wfn1_wave,
     &                               en1a,en1b,ec1a,ec1b)
      implicit none
C>
C> \brief Print the wavefunction analysis
C>
#include "stdio.fh"
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
!>    The natural orbital occupation numbers
      double precision, allocatable   :: occa(:), occb(:)
      double precision, allocatable   :: ci(:), cj(:)
c
      double precision :: tanalyze
      double precision :: valmin
      integer          :: imin
c
      integer, allocatable :: irrep(:)
c
      integer :: ii, jj   ! counters
      integer :: ilo, ihi ! range limits
      integer :: nmo      ! local copy of number of orbitals
c
      character*255 blob
c
      tanalyze = 0.025d0
      ilo = 1
      ihi = wfn1_param%nmo
      nmo = wfn1_param%nmo
c
c     allocate(e1a(1:wfn1_param%nmo),e1b(1:wfn1_param%nmo))
      allocate(occa(1:nmo),occb(1:nmo))
      allocate(ci(1:nmo),cj(1:nmo))
      allocate(irrep(1:nmo))
c     e1a   = 0.0d0
c     e1b   = 0.0d0
      irrep = 0
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,wfn1_param%nea,
     &                   wfn1_wave%icnfa,wfn1_wave%dcnta,occa)
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,wfn1_param%neb,
     &                   wfn1_wave%icnfb,wfn1_wave%dcntb,occb)
c
c     call wfn1_print_occ(nmo,occa,occb)
c
      blob = "WFN1 Final Natural orbital energies and occupations"
      write(LuOut,*)
      call util_print_centered(LuOut,blob, 40, .true.)
      write(LuOut,*)
      write(LuOut,'(12x,"  Alpha ",4x,"  Beta  ",8x,
     &                  "  Alpha ",4x,"  Beta  ")')
      write(LuOut,'(12x,"--------",4x,"--------",8x,
     &                  "--------",4x,"--------")')
      do ii = 1, nmo
        write(LuOut,'(i8,f12.4,f12.4,4x,f12.4,f12.4)')ii,
     &        en1a(ii),en1b(ii),occa(ii),occb(ii)
     &
      enddo
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
c
c     Beta
c
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
      call ga_sync()
c
      blob = "WFN1 Final Correlation function energies and occupations"
      write(LuOut,*)
      call util_print_centered(LuOut,blob, 40, .true.)
      write(LuOut,*)
      write(LuOut,'(12x,"  Alpha ",4x,"  Beta  ",8x,
     &                  "  Alpha ",4x,"  Beta  ")')
      write(LuOut,'(12x,"--------",4x,"--------",8x,
     &                  "--------",4x,"--------")')
      do ii = ilo, ihi
        write(LuOut,'(i8,f12.4,f12.4,4x,f12.4,f12.4)')ii,
     &        ec1a(ii),ec1b(ii),occa(ii),occb(ii)
     &
      enddo
c
      blob = "WFN1 Final alpha Correlation Function Analysis"
      write(LuOut,*)
      call util_print_centered(LuOut,blob, 40, .true.)
      write(LuOut,*)
      call util_flush(LuOut)
      call ga_print(wfn1_wave%g_ca)
      call util_flush(LuOut)
c
      blob = "WFN1 Final beta Correlation Function Analysis"
      write(LuOut,*)
      call util_print_centered(LuOut,blob, 40, .true.)
      write(LuOut,*)
      call util_flush(LuOut)
      call ga_print(wfn1_wave%g_cb)
      call util_flush(LuOut)
c
      blob = "WFN1 Final Alpha Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_na,blob,.true.,ec1a,.false.,
     &                       irrep,.true.,occa)
c
      blob = "WFN1 Final Beta Natural Orbital Analysis"
      call movecs_print_anal(wfn1_inst%basis,ilo,ihi,tanalyze,
     &                       wfn1_wave%g_nb,blob,.true.,ec1b,.false.,
     &                       irrep,.true.,occb)
c
      deallocate(occa,occb)
      deallocate(ci,cj)
c     deallocate(e1a,e1b)
c
      end
c
c-----------------------------------------------------------------------
