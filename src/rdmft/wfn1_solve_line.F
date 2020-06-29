c-----------------------------------------------------------------------
c
C> \brief Solve for the minimal energy using line searches
C>
C> The Fock matrices are used as search directions and the energy is
C> minimized along these directions. In practice a quadratic expression
C> is fitted to the energy at different step sizes and the step size with
C> minimum energy is ultimately picked as the solution. 
C> Multiple line searches are issued until overal convergence is
C> reached.
C>
      logical
     &function wfn1_solve_line(wfn1_param,wfn1_inst,wfn1_mat,
     &                         wfn1_wave,wfn1_deriv,wfn1_energy,
     &                         wfn1_trial_wv,wfn1_trial_drv,
     &                         wfn1_part_drv,
     &                         g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &                         g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                         v_dna,v_dnb,v_dca,v_dcb,
     &                         v_dnna,v_dnnb,v_dcna,v_dcnb,
     &                         oprerr,curit)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_size.fh"
#include "wfn1_trace.fh"
#include "wfn1_en.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
c
c     Arguments
c
!     Parameter instances (const complex data structures)
      type(wfn1_prminst), intent(in) :: wfn1_inst
!     Calculation parameters
      type(wfn1_prm), intent(in)     :: wfn1_param
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
      integer :: g_bfbf
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
!     The current iteration count
      integer :: curit
      integer :: curit0
c
      double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnnb(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcnb(1:wfn1_param%nmo)
c
c     Local
c
      integer :: iproc
      integer :: maxit
      logical :: oconverged
c
      double precision :: conv_tol 
c
      type(wfn1_wfn)   :: wfn1_t_wv
      type(wfn1_drv)   :: wfn1_t_drv
      type(wfn1_error) :: error
c
      double precision :: error_nst ! Natural orbital error at start
      double precision :: error_cst ! Correlation function error at start
c
      integer, parameter :: DIAG_N = 1 ! diagonalize nat orb Fock
      integer, parameter :: DIAG_C = 2 ! diagonalize cor fnc Fock
      integer  :: diag_NorC = DIAG_C
c
      iproc    = ga_nodeid()
      conv_tol = wfn1_param%convergence
      maxit    = wfn1_param%maxiter
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
     &              g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &              g_bfmo,
     &              g_momo1,g_momo2,g_momo3,g_momo4,
     &              v_dna,v_dnb,v_dca,v_dcb,
     &              v_dnna,v_dnnb,v_dcna,v_dcnb)
c     call wfn1_skew(wfn1_param,wfn1_deriv,v_dna,v_dnb,v_dca,v_dcb)
c     call wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_deriv)
      error      = wfn1_trace(wfn1_deriv,wfn1_deriv)
      oconverged = (error%total < conv_tol)
      curit      = 0
      curit0     = 0
      error_nst  = min(error%errna+error%errnb,error%errca+error%errcb)
      error_cst  = min(error%errna+error%errnb,error%errca+error%errcb)
      if (iproc.eq.0) then
        if (.not.oprerr) then
          write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                    error%total
        else
          write(luout,'(i10,f20.8,5f16.6)')curit,wfn1_energy%E_tot,
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
        if      (error%errca+error%errcb > 0.5d0*conv_tol) then
          diag_NorC = diag_C
        else if (error%errna+error%errnb > 0.5d0*conv_tol) then
          diag_NorC = diag_N
        endif
        if      (diag_NorC.eq.diag_N) then
          error_cst = error%errca+error%errcb
          call ga_zero(wfn1_deriv%g_dca)
          call ga_zero(wfn1_deriv%g_dcb)
        else if (diag_NorC.eq.diag_C) then
          error_nst = error%errna+error%errnb
          call ga_zero(wfn1_deriv%g_dna)
          call ga_zero(wfn1_deriv%g_dnb)
        endif
c       if (error%errna+error%errnb > 0.5d0*conv_tol .and.
c    &      error%errca+error%errcb > 0.5d0*conv_tol) then
c         if      (diag_NorC.eq.diag_N) then
c           if (error%errna+error%errnb < 0.1d0*error_nst .or.
c    &          curit-curit0 > 2) then
c             diag_NorC = diag_C
c             curit0    = curit
c           endif
c         else if (diag_NorC.eq.diag_C) then
c           if (error%errca+error%errcb < 0.1d0*error_cst .or.
c    &          curit-curit0 > 2) then
c             diag_NorC = diag_N
c             curit0    = curit
c           endif
c         endif
c         if      (diag_NorC.eq.diag_N) then
c           error_cst = error%errca+error%errcb
c           call ga_zero(wfn1_deriv%g_dca)
c           call ga_zero(wfn1_deriv%g_dcb)
c         else if (diag_NorC.eq.diag_C) then
c           error_nst = error%errna+error%errnb
c           call ga_zero(wfn1_deriv%g_dna)
c           call ga_zero(wfn1_deriv%g_dnb)
c         endif
c       else if (error%errna+error%errnb >
c    &           error%errca+error%errcb) then
c         diag_NorC = diag_N
c         call ga_zero(wfn1_deriv%g_dca)
c         call ga_zero(wfn1_deriv%g_dcb)
c       else
c         diag_NorC = diag_C
c         call ga_zero(wfn1_deriv%g_dna)
c         call ga_zero(wfn1_deriv%g_dnb)
c       endif
c
c       end of alternating optimization
c
        call wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       wfn1_deriv,wfn1_energy,wfn1_trial_wv,
     &                       wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &                       g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &                       g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                       v_dna,v_dnb,v_dca,v_dcb,
     &                       v_dnna,v_dnnb,v_dcna,v_dcnb)
c
        error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
        oconverged = (error%total < conv_tol)
c
c       swap trial wavefunctions, derivatives and energies with the 
c       regular ones.
c
        wfn1_t_wv      = wfn1_trial_wv
        wfn1_trial_wv  = wfn1_wave
        wfn1_wave      = wfn1_t_wv
c
        wfn1_t_drv     = wfn1_trial_drv
        wfn1_trial_drv = wfn1_deriv
        wfn1_deriv     = wfn1_t_drv
c
        wfn1_energy    = wfn1_trial_en

        if (iproc.eq.0) then
          if (.not.oprerr) then
            write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                      error%total
          else
            if      (diag_NorC.eq.diag_C) then
              write(luout,'(i10,f20.8,5f16.6,a2)')
     &              curit,wfn1_energy%E_tot,
     &              error%total,
     &              error%errna,error%errnb,
     &              error%errca,error%errcb,'C'
            else if (diag_NorC.eq.diag_N) then
              write(luout,'(i10,f20.8,5f16.6,a2)')
     &              curit,wfn1_energy%E_tot,
     &              error%total,
     &              error%errna,error%errnb,
     &              error%errca,error%errcb,'N'
            else
              call errquit("should never get here",0,UERR)
            endif
          endif
        endif
      enddo
      wfn1_solve_line = oconverged
      end function wfn1_solve_line
c
c-----------------------------------------------------------------------
