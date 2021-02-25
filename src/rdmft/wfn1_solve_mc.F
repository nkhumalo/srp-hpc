c-----------------------------------------------------------------------
c
C> \brief Solve for the minimal energy using Monte-Carlo updates
C>
C> Use Monte-Carlo to minimize the energy. This is obviously a dumb
C> approach as it is very expensive. However, this approach does not
C> require proper Fock matrices to work. So it should provide a good 
C> test (albeit slow).
C>
      logical
     &function wfn1_solve_mc(wfn1_param,wfn1_inst,wfn1_mat,
     &                       wfn1_wave,
!    &                       wfn1_deriv,
     &                       wfn1_energy,wfn1_trial_wv,
!    &                       wfn1_trial_drv,wfn1_part_drv,
     &                       g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &                       g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                       g_4t1,g_4t2,
     &                       gz_momo1,
!    &                       v_dna,v_dnb,v_dca,v_dcb,
!    &                       v_dnna,v_dnnb,v_dcna,v_dcnb,
     &                       oprerr,curit)
      implicit none
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_size.fh"
!include "wfn1_trace.fh"
!include "wfn1_en.fh"
#include "errquit.fh"
#include "global.fh"
#include "stdio.fh"
c
c     Arguments
c
!     Parameter instances (const complex data structures)
      type(wfn1_prminst), intent(in) :: wfn1_inst
!     Calculation parameters
      type(wfn1_prm), intent(inout)  :: wfn1_param
!     Matrix constants
      type(wfn1_matcnst), intent(in) :: wfn1_mat
!     The current wave function
      type(wfn1_wfn), intent(inout)  :: wfn1_wave
!     The current wave function derivative
!     type(wfn1_drv), intent(inout)  :: wfn1_deriv
!     The energy terms
      type(wfn1_en), intent(inout)   :: wfn1_energy
c
!     The trial wave function
      type(wfn1_wfn), intent(inout)  :: wfn1_trial_wv
!     The trial wave function derivative
!     type(wfn1_drv), intent(inout)  :: wfn1_trial_drv
!     The non-correlation, non-thermal part of wave function derivative
!     type(wfn1_drv), intent(inout)  :: wfn1_part_drv
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
!     GA dimensions nbf x nbf x nbf x nbf
      integer :: g_4t1
!     GA dimensions nbf x nbf x nbf x nbf
      integer :: g_4t2
!     GA dimensions nmo x nmo complex
      integer :: gz_momo1
c
!     Print detailed errors or not
      logical :: oprerr
!     The current iteration count
      integer :: curit
c
!     double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dnna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dnnb(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcnb(1:wfn1_param%nmo)
c
c     Local
c
      integer :: nmo
      integer :: iproc
      integer :: maxit
c
      double precision :: rlimo ! the old rlimit
c
      type(wfn1_wfn)   :: wfn1_t_wv
      type(wfn1_error) :: error
c
      iproc = ga_nodeid()
      maxit = wfn1_param%maxiter
      rlimo = wfn1_param%rlim
      nmo   = wfn1_param%nmo
      wfn1_solve_mc = .false.
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
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &              wfn1_deriv,wfn1_part_drv,
     &              wfn1_energy,
     &              g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &              g_bfmo,
     &              g_momo1,g_momo2,g_momo3,g_momo4,
     &              g_4t1,g_4t2)
!    &              v_dna,v_dnb,v_dca,v_dcb,
!    &              v_dnna,v_dnnb,v_dcna,v_dcnb)
!     error       = wfn1_trace(wfn1_deriv,wfn1_deriv)
!
!     Estimate of the distance to the solution based on the hypercube
!     defined by the interval in which the random numbers are chosen.
!     The estimate corresponds to the distance of the center to one
!     of the corners of the hypercube.
!
      error%total = 2.0d0*rlimo*nmo
      curit       = 0
      if (iproc.eq.0) then
        if (.not.oprerr) then
          write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                    error%total
        else
          write(luout,'(i10,f20.8,5f16.6)')curit,wfn1_energy%E_tot,
     &                                     error%total,
     &                                     error%errna,error%errnb,
     &                                     abs(error%errca),
     &                                     abs(error%errcb)
        endif
      endif
      do while (curit.lt.maxit.and.
     &          error%total.gt.wfn1_param%convergence)
        curit = curit + 1
        call wfn1_mcsearch(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                     wfn1_deriv,
     &                     wfn1_energy,wfn1_trial_wv,
!    &                     wfn1_trial_drv,wfn1_part_drv,
     &                     wfn1_trial_en,
     &                     g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &                     g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                     g_4t1,g_4t2,
     &                     gz_momo1)
!    &                     v_dna,v_dnb,v_dca,v_dcb,
!    &                     v_dnna,v_dnnb,v_dcna,v_dcnb)
c
!       error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
        if (wfn1_param%rlim.lt.rlimo) then
          error%total = 2.0d0*rlimo*nmo
          rlimo       = wfn1_param%rlim
        endif
c
c       swap trial wavefunctions, derivatives and energies with the 
c       regular ones.
c
        wfn1_t_wv      = wfn1_trial_wv
        wfn1_trial_wv  = wfn1_wave
        wfn1_wave      = wfn1_t_wv
c
        wfn1_energy    = wfn1_trial_en
c
        if (iproc.eq.0) then
          if (.not.oprerr) then
            write(luout,'(i10,f20.8,f16.6)')curit,wfn1_energy%E_tot,
     &                                      error%total
          else
            write(luout,'(i10,f20.8,5f16.6)')
     &            curit,wfn1_energy%E_tot,
     &            error%total,
     &            error%errna,error%errnb,
     &            abs(error%errca),abs(error%errcb)
          endif
        endif
      enddo
      wfn1_solve_mc = .true.
c
      end
c
c-----------------------------------------------------------------------
