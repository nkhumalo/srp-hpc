c-----------------------------------------------------------------------
c
!> \brief Store the wave function
!>
!> Store the wave function using the conventional I/O routines of 
!> NWChem for this purpose. Other code modules in NWChem will not be
!> aware of the correlation functions and hence they will be stored 
!> on a separate file. The natural orbitals stored this way should be
!> available for other modules to use.
!>
!> This routine is essentially implemented as a wrapper around the
!> conventional `movecs_write` routine.
!>
      subroutine wfn1_write_movecs(wfn1_param,wfn1_inst,wfn1_wave,
     &                             v_dna,v_dnb,v_dca,v_dcb)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
!>    The parameters of the calculation
      type(wfn1_prm),     intent(in) :: wfn1_param
!>    The parameter instances of the calculation
      type(wfn1_prminst), intent(in) :: wfn1_inst
!>    The wave function
      type(wfn1_wfn),     intent(in) :: wfn1_wave
!>    The \f$\alpha\f$ natural orbital energies
      double precision,   intent(in) :: v_dna(wfn1_param%nmo)
!>    The \f$\beta\f$ natural orbital energies
      double precision,   intent(in) :: v_dnb(wfn1_param%nmo)
!>    The \f$\alpha\f$ correlation function energies
      double precision,   intent(in) :: v_dca(wfn1_param%nmo)
!>    The \f$\beta\f$ correlation function energies
      double precision,   intent(in) :: v_dcb(wfn1_param%nmo)
c
      double precision, allocatable :: ocn(:,:)
      double precision, allocatable :: occ(:,:)
      double precision, allocatable :: evn(:,:)
      double precision, allocatable :: evc(:,:)
c
      integer :: g_vecs(2)
      integer :: nmo(2)
      integer :: ii        ! counter
c
      character*256 covecs_out
c
      allocate(ocn(wfn1_param%nmo,2))
      allocate(occ(wfn1_param%nmo,2))
      allocate(evn(wfn1_param%nmo,2))
      allocate(evc(wfn1_param%nmo,2))
      nmo(1) = wfn1_param%nmo
      nmo(2) = wfn1_param%nmo
c
c     Set the occupation numbers up
c
      call wfn1_calc_occ(wfn1_wave%g_ca,wfn1_param%nmo,wfn1_param%nea,
     &                   wfn1_wave%icnfa,wfn1_wave%dcnta,ocn(:,1))
      call wfn1_calc_occ(wfn1_wave%g_cb,wfn1_param%nmo,wfn1_param%neb,
     &                   wfn1_wave%icnfb,wfn1_wave%dcntb,ocn(:,2))
      occ = 0.0d0
      do ii = 1, wfn1_param%nea
        occ(wfn1_wave%icnfa(ii),1) = wfn1_wave%dcnta(ii)
      enddo
      do ii = 1, wfn1_param%neb
        occ(wfn1_wave%icnfb(ii),2) = wfn1_wave%dcntb(ii)
      enddo
c
c     Set the energies up
c
      evn(:,1) = v_dna(:)
      evn(:,2) = v_dnb(:)
      evc(:,1) = v_dca(:)
      evc(:,2) = v_dcb(:)
c
c     Store the correlation functions
c
      g_vecs(1)  = wfn1_wave%g_ca
      g_vecs(2)  = wfn1_wave%g_cb
      covecs_out = wfn1_param%movecs_out
      ii = len(trim(covecs_out))
      if (covecs_out(ii-5:ii)=="movecs") then
        covecs_out = covecs_out(1:ii-6)//"covecs"
      else
        covecs_out = covecs_out(1:ii)//".covecs"
      endif
      call movecs_write(wfn1_inst%rtdb,wfn1_inst%basis,
     &                  trim(covecs_out),
     &                  "WFN1","",wfn1_param%nmo,2,nmo,
     &                  occ,wfn1_param%nmo,evc,wfn1_param%nmo,g_vecs)
c
c     Store the natural orbitals
c
      g_vecs(1) = wfn1_wave%g_na
      g_vecs(2) = wfn1_wave%g_nb
      call movecs_write(wfn1_inst%rtdb,wfn1_inst%basis,
     &                  trim(wfn1_param%movecs_out),
     &                  "WFN1","",wfn1_param%nmo,2,nmo,
     &                  ocn,wfn1_param%nmo,evn,wfn1_param%nmo,g_vecs)
c
      deallocate(evn)
      deallocate(evc)
      deallocate(occ)
      deallocate(ocn)
c
      end subroutine wfn1_write_movecs
c
c-----------------------------------------------------------------------
