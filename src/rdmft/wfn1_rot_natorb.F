c-----------------------------------------------------------------------
c
C> \brief Rotate natural orbitals
C>
C> Update the natural orbitals in the context of a finite temperature
C> calculation. Stop when either convergence is reached or a small 
C> number of iterations are completed.
C>
C> The aim is to move the natural orbitals away from orbitals 
C> corresponding to an idem-potent 1-electron density matrix.
C>
      subroutine wfn1_rot_natorb(wfn1_param_in,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,wfn1_deriv,wfn1_energy,
     &                           wfn1_trial_wv,wfn1_trial_drv,
     &                           wfn1_part_drv,
     &                           g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                           g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                           v_dna,v_dnb,v_dca,v_dcb,
     &                           v_dnna,v_dnnb,v_dcna,v_dcnb,
     &                           oprerr,curit)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_size.fh"
#include "wfn1_trace.fh"
#include "wfn1_en.fh"
#include "global.fh"
#include "errquit.fh"
#include "stdio.fh"
c
c     Arguments
c
!     Parameter instances (const complex data structures)
      type(wfn1_prminst), intent(in) :: wfn1_inst
!     Calculation parameters
      type(wfn1_prm), intent(in)     :: wfn1_param_in
!     Matrix constants
      type(wfn1_matcnst), intent(in) :: wfn1_mat
!     The current wave function
      type(wfn1_wfn), intent(inout)  :: wfn1_wave
!     The current wave function derivative
      type(wfn1_drv), intent(inout)  :: wfn1_deriv
!     The energy terms
      type(wfn1_en), intent(inout)   :: wfn1_energy
c
!     The trial wave function
      type(wfn1_wfn), intent(inout)  :: wfn1_trial_wv
!     The trial wave function derivative
      type(wfn1_drv), intent(inout)  :: wfn1_trial_drv
!     The non-correlation, non-thermal part of wave function derivative
      type(wfn1_drv), intent(inout)  :: wfn1_part_drv
!     The trial energy terms
      type(wfn1_en)      :: wfn1_trial_en  ! the trial energy
!     GA dimensions nbf x nbf
      integer :: g_bfbf1
!     GA dimensions nbf x nbf
      integer :: g_bfbf2
!     GA dimensions nbf x nbf
      integer :: g_bfbf3
!     GA dimensions nbf x nbf
      integer :: g_bfbf4
!     GA dimensions nbf x nmo
      integer :: g_bfmo
!     GA dimensions nmo x nmo
      integer :: g_momo1
!     GA dimensions nmo x nmo
      integer :: g_momo2
!     GA dimensions nmo x nmo
      integer :: g_momo3
!     GA dimensions nmo x nmo
      integer :: g_momo4
c
!     Print detailed errors or not
      logical :: oprerr
c
      double precision, intent(out) :: v_dna(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dnb(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dca(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dcb(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dnna(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dnnb(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dcna(1:wfn1_param_in%nmo)
      double precision, intent(out) :: v_dcnb(1:wfn1_param_in%nmo)
c
c     Local
c
      integer :: iproc
      integer :: maxit
      integer :: curit
      integer :: nmo
      integer :: nbf
      integer :: ii ! counter
      logical :: oconverged
c
      integer, allocatable :: icola(:)
      integer, allocatable :: icolb(:)
c
      double precision :: conv_tol 
c
      type(wfn1_wfn)   :: wfn1_t_wv
      type(wfn1_drv)   :: wfn1_t_drv
      type(wfn1_prm)   :: wfn1_param
      type(wfn1_error) :: error
c
      integer, parameter :: DIAG_N = 1 ! diagonalize nat orb Fock
      integer, parameter :: DIAG_C = 2 ! diagonalize cor fnc Fock
      integer  :: diag_NorC = DIAG_N
cDEBUG
c     double precision :: ea, eb
cDEBUG
c
      wfn1_param           = wfn1_param_in
      wfn1_param%maxiter   = 20
      wfn1_param%Tbath     = wfn1_param%Tinitial
      wfn1_param%corr_expr = "wfn1_mx"
c
      iproc      = ga_nodeid()
      conv_tol   = wfn1_param%convergence
      maxit      = wfn1_param%maxiter
      nmo        = wfn1_param%nmo
      nbf        = wfn1_param%nbf
c
      allocate(icola(1:nmo))
      allocate(icolb(1:nmo))
c
      if (iproc.eq.0) then
        if (.not.oprerr) then
          write(luout,'(a10,a20,a16)')"iteration","energy","error"
          write(luout,'(a10,a20,a16)')"---------","------","-----"
        else
          write(luout,'(a10,a20,a16,4a16)')"iteration","energy","error",
     &          "err Na","err Nb","err Ca","err Cb"
          write(luout,'(a10,a20,a16,4a16)')"---------","------","-----",
     &          "------","------","------","------"
        endif
      endif
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,wfn1_deriv,
     &              wfn1_part_drv,wfn1_energy,
     &              g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &              g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &              v_dna,v_dnb,v_dca,v_dcb,
     &              v_dnna,v_dnnb,v_dcna,v_dcnb)
      call wfn1_symm(wfn1_param,wfn1_deriv,v_dna,v_dnb,v_dca,v_dcb)
      error      = wfn1_trace(wfn1_deriv,wfn1_deriv)
      oconverged = (error%total < conv_tol)
      curit      = 0
      if (iproc.eq.0) then
        if (.not.oprerr) then
          write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                    error%total
        else
          write(luout,'(i10,f20.8,5e16.6)')curit,wfn1_energy%E_tot,
     &                                     error%total,
     &                                     error%errna,error%errnb,
     &                                     error%errca,error%errcb
        endif
      endif
      do while ((.not.oconverged).and.(curit.lt.maxit))
        curit = curit + 1
c
c       Alternating optimization of natural orbitals and correlation
c       functions for better convergence.
c
c       if (error%errna+error%errnb > 0.5d0*conv_tol .and.
c    &      error%errca+error%errcb > 0.5d0*conv_tol) then
c         if (mod(curit,2).eq.0) then
c           diag_NorC = DIAG_N
c         else
c           diag_NorC = DIAG_C
c         endif
c       else if (error%errna+error%errnb >
c    &           error%errca+error%errcb) then
c         diag_NorC = DIAG_N
c       else
c         diag_NorC = DIAG_C
c       endif
c
c       Different convergence scheme where the natural orbitals are
c       converged first, and only then the correlation functions
c       are converged.
c
c       if (error%errna+error%errnb > 0.5d0*conv_tol) then
c         diag_NorC = DIAG_N
c       else
c         diag_NorC = DIAG_C
c       endif
c
c       Different convergence scheme where the correlation functions are
c       converged first, and only then the natural orbitals
c       are converged.
c
        diag_NorC = DIAG_N
c
c       end of alternating optimization
c
        if      (diag_NorC.eq.DIAG_C) then
cDEBUG
c         do ii = 1, nmo
c           call ga_get(wfn1_deriv%g_dca,ii,ii,ii,ii,ea,1)
c           call ga_get(wfn1_deriv%g_dcb,ii,ii,ii,ii,eb,1)
c           write(*,'(" HVD: eorb:",i4,2e18.6)')ii,ea,eb
c         enddo
cDEBUG
          call wfn1_scale_offdiag(wfn1_deriv%g_dca,wfn1_param%damp_c)
          call wfn1_scale_offdiag(wfn1_deriv%g_dcb,wfn1_param%damp_c)
          call ga_diag_std(wfn1_deriv%g_dca,g_momo1,v_dca)
          call ga_diag_std(wfn1_deriv%g_dcb,g_momo2,v_dcb)
          call wfn1_column_max(wfn1_param,g_momo1,icola)
          call wfn1_column_max(wfn1_param,g_momo2,icolb)
          call wfn1_backorder(wfn1_param,g_momo1,icola,g_momo3)
          call wfn1_backorder(wfn1_param,g_momo2,icolb,g_momo4)
c         call ga_copy(g_momo1,g_momo3)
c         call ga_copy(g_momo2,g_momo4)
          call ga_dgemm('N','N',nmo,nmo,nmo,
     &                  1.0d0,wfn1_wave%g_ca,g_momo3,0.0d0,g_momo1)
          call ga_dgemm('N','N',nmo,nmo,nmo,
     &                  1.0d0,wfn1_wave%g_cb,g_momo4,0.0d0,g_momo2)
          call ga_copy(g_momo1,wfn1_wave%g_ca)
          call ga_copy(g_momo2,wfn1_wave%g_cb)
        else if (diag_NorC.eq.DIAG_N) then
          call wfn1_scale_offdiag(wfn1_deriv%g_dna,wfn1_param%damp_n)
          call wfn1_scale_offdiag(wfn1_deriv%g_dnb,wfn1_param%damp_n)
          call ga_diag_std(wfn1_deriv%g_dna,g_momo1,v_dna)
          call ga_diag_std(wfn1_deriv%g_dnb,g_momo2,v_dnb)
          call wfn1_column_max(wfn1_param,g_momo1,icola)
          call wfn1_column_max(wfn1_param,g_momo2,icolb)
          call wfn1_backorder(wfn1_param,g_momo1,icola,g_momo3)
          call wfn1_backorder(wfn1_param,g_momo2,icolb,g_momo4)
          call ga_dgemm('N','N',nmo,nmo,nmo,
     &                  1.0d0,wfn1_wave%g_na,g_momo3,0.0d0,g_momo1)
          call ga_dgemm('N','N',nmo,nmo,nmo,
     &                  1.0d0,wfn1_wave%g_nb,g_momo4,0.0d0,g_momo2)
          call ga_copy(g_momo1,wfn1_wave%g_na)
          call ga_copy(g_momo2,wfn1_wave%g_nb)
        else
          call errquit("wfn1_solve_diag: invalid diagonalization",
     &                 diag_NorC,UERR)
        endif
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                wfn1_deriv,wfn1_part_drv,wfn1_energy,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,
     &                g_momo1,g_momo2,g_momo3,g_momo4,
     &                v_dna,v_dnb,v_dca,v_dcb,
     &                v_dnna,v_dnnb,v_dcna,v_dcnb)
        call wfn1_symm(wfn1_param,wfn1_deriv,v_dna,v_dnb,v_dca,v_dcb)
        error = wfn1_trace(wfn1_deriv,wfn1_deriv)
        oconverged = (error%total < conv_tol)
c
        if (iproc.eq.0) then
          if (.not.oprerr) then
            write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                      error%total
          else
            if      (diag_NorC.eq.diag_C) then
              write(luout,'(i10,f20.8,5e16.6,a2)')
     &              curit,wfn1_energy%E_tot,error%total,
     &              error%errna,error%errnb,error%errca,error%errcb,'C'
            else if (diag_NorC.eq.diag_N) then
              write(luout,'(i10,f20.8,5e16.6,a2)')
     &              curit,wfn1_energy%E_tot,error%total,
     &              error%errna,error%errnb,error%errca,error%errcb,'N'
            else
              call errquit("should never get here",0,UERR)
            endif
          endif
        endif
      enddo
      deallocate(icolb)
      deallocate(icola)
      end subroutine wfn1_rot_natorb
c
c-----------------------------------------------------------------------
