c-----------------------------------------------------------------------
c
!> \brief Calculate the redistribution energy
!>
!> Correlation effects change the distribution of electrons across the
!> orbitals. As the electrons typically occupy the lowest energy 
!> orbitals redistributing the electrons  costs some energy. This 
!> subroutine calculates this energy costs as well as the derivatives
!> of this energy.
!>
!> The redistribution energy is simply calculated by first calculating
!> the energy of the reference state, then calculating the energy of
!> the current state and subsequently the difference. 
!> The energy of any particular state is calculated as
!> \f{eqnarray*}{
!>    E_{ref} &=& \sum_i \epsilon_i \rho^{ref}_i \\\\
!>    E_{state} &=& \sum_i \epsilon_i \rho_i \\\\
!>    E_{redist} &=& E_{state} - E_{ref}
!> \f}
!> where \f$\rho^{ref}_i\f$ are simply the occupation numbers of the
!> correlation free reference state. I.e. they are either \f$0\f$ 
!> or \f$1\f$ depending on the state.
!> The \f$\rho_i\f$ are the actual occupation numbers of the correlated
!> state. 
!>
!> Hence we have that if there is no correlation the redistribution
!> energy is \f$0\f$. Also note that differentating this term only
!> involves differentiating \f$E_{state}\f$ as the occupation numbers
!> in \f$E_{ref}\f$ are simply constants.
!>
!> The current implementation is valid for ground states only.
!>
      subroutine wfn1_e_g_redist(wfn1_param,wfn1_wave,v_dca,v_dcb,
     &                           E_redist)
      use wfn1_fac
      implicit none
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
c
      type(wfn1_prm),   intent(in)    :: wfn1_param
      type(wfn1_wfn),   intent(in)    :: wfn1_wave
c
      double precision, intent(inout) :: E_redist
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer :: nmo
      integer :: nea
      integer :: neb
      integer :: ii
      double precision, allocatable :: ea(:)
      double precision, allocatable :: eb(:)
      double precision, allocatable :: occa(:)
      double precision, allocatable :: occb(:)
c
!     double precision, external :: wfn1_redist
c
      if (wfn1_redist().le.0.0d0) return
c
      nmo = wfn1_param%nmo
      nea = wfn1_param%nea
      neb = wfn1_param%neb
c
c     Account for the reference energy
c
      allocate(ea(1:nmo))
      allocate(eb(1:nmo))
      ea = orb_ea
      eb = orb_eb
      call wfn1_sort(ea,nmo)
      call wfn1_sort(eb,nmo)
      do ii = 1, nea
        E_redist = E_redist - wfn1_redist()*ea(ii)
      enddo
      do ii = 1, neb
        E_redist = E_redist - wfn1_redist()*eb(ii)
      enddo
      deallocate(eb)
      deallocate(ea)
c
c     Calculate the present state energy
c
      allocate(occa(1:nmo))
      allocate(occb(1:nmo))
      occa = 0.0d0
      occb = 0.0d0
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
      do ii = 1, nmo
        E_redist = E_redist + wfn1_redist()*orb_ea(ii)*occa(ii)
      enddo
      do ii = 1, nmo
        E_redist = E_redist + wfn1_redist()*orb_eb(ii)*occb(ii)
      enddo
      deallocate(occb)
      deallocate(occa)
c
c     Add the terms to the Fock matrix
c
      v_dca = v_dca + wfn1_redist()*orb_ea
      v_dcb = v_dcb + wfn1_redist()*orb_eb
c
      end subroutine wfn1_e_g_redist
c
c-----------------------------------------------------------------------
