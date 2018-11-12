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
     &                         g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                         v_dna,v_dnb,v_dca,v_dcb,
     &                         oprerr,curit)
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_size.fh"
#include "wfn1_trace.fh"
#include "wfn1_en.fh"
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
!     The trial energy terms
      type(wfn1_en)      :: wfn1_trial_en  ! the trial energy
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
c
      double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
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
     &              wfn1_energy,g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &              v_dna,v_dnb,v_dca,v_dcb)
      call wfn1_skew(wfn1_param,wfn1_deriv,v_dna,v_dnb,v_dca,v_dcb)
c     call wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_deriv)
      error      = wfn1_trace(wfn1_deriv,wfn1_deriv)
      oconverged = (error%total < conv_tol)
      curit      = 0
      do while ((.not.oconverged).and.(curit.lt.maxit))
        curit = curit + 1
c
c       Alternating optimization of natural orbitals and correlation
c       functions for better convergence.
c
        if (error%errna+error%errnb > 0.5d0*conv_tol .and.
     &      error%errca+error%errcb > 0.5d0*conv_tol) then
          if (mod(curit,2).eq.0) then
            call ga_zero(wfn1_deriv%g_dca)
            call ga_zero(wfn1_deriv%g_dcb)
          else
            call ga_zero(wfn1_deriv%g_dna)
            call ga_zero(wfn1_deriv%g_dnb)
          endif
        else if (error%errna+error%errnb >
     &           error%errca+error%errcb) then
          call ga_zero(wfn1_deriv%g_dca)
          call ga_zero(wfn1_deriv%g_dcb)
        else
          call ga_zero(wfn1_deriv%g_dna)
          call ga_zero(wfn1_deriv%g_dnb)
        endif
c
c       end of alternating optimization
c
        call wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       wfn1_deriv,wfn1_energy,
     &                       wfn1_trial_wv,wfn1_trial_drv,wfn1_trial_en,
     &                       g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                       v_dna,v_dnb,v_dca,v_dcb)
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
            write(luout,'(i10,f20.8,5f16.6)')curit,wfn1_energy%E_tot,
     &                                       error%total,
     &                                       error%errna,error%errnb,
     &                                       error%errca,error%errcb
          endif
        endif
      enddo
      wfn1_solve_line = oconverged
      end function wfn1_solve_line
c
c-----------------------------------------------------------------------
