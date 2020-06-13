c-----------------------------------------------------------------------
c
C>
C> \brief Perform a linesearch for the minimum total energy
C>
C> Given the gradient of the energy wrt. the wave function parameters
C> this gradient is used as the search direction to find the minimum
C> of the total energy. The basic approach is to start from the
C> energy and gradient at the initial, and the energy after taking a
C> step, and fit a parabola through these three data points. Looking
C> for the minimum of that parabola a suggestion for the minimum is
C> obtained. In subsequent steps the energy at three points might be
C> used to fit a parabola. 
C>
C> When enough progress has been made the subroutine returns the 
C> current energies, the current wavefunction, and the current gradients
C> for the next linesearch (or the final solution if the calculation
C> is converged).
C>
      subroutine wfn1_linesearch(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,wfn1_deriv,wfn1_energy,
     &                           wfn1_trial_wv,wfn1_trial_drv,
     &                           wfn1_part_drv,wfn1_trial_en,
     &                           g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                           g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                           v_dna,v_dnb,v_dca,v_dcb,
     &                           v_dnna,v_dnnb,v_dcna,v_dcnb)
      implicit none
#include "global.fh"
#include "errquit.fh"
#include "wfn1_trace.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_drv.fh"
#include "wfn1_en.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(inout) :: wfn1_wave
      type(wfn1_drv),     intent(inout) :: wfn1_deriv
      type(wfn1_en),      intent(inout) :: wfn1_energy
      type(wfn1_wfn),     intent(inout) :: wfn1_trial_wv
      type(wfn1_drv),     intent(inout) :: wfn1_trial_drv
      type(wfn1_drv),     intent(inout) :: wfn1_part_drv
      type(wfn1_en),      intent(inout) :: wfn1_trial_en
c
      integer, intent(inout) :: g_bfbf1 ! nbf x nbf scratch GA
      integer, intent(inout) :: g_bfbf2 ! nbf x nbf scratch GA
      integer, intent(inout) :: g_bfbf3 ! nbf x nbf scratch GA
      integer, intent(inout) :: g_bfbf4 ! nbf x nbf scratch GA
      integer, intent(inout) :: g_bfmo  ! nbf x nmo scratch GA
      integer, intent(inout) :: g_momo1 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_momo2 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_momo3 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_momo4 ! nmo x nmo scratch GA
c
c     Natural orbital energies
c
      double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
c
c     Correlation function energies
c
      double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
c
c     Natural orbital expectation values of the Correlation
c     function Fock matrices
c
      double precision, intent(out) :: v_dnna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dnnb(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcna(1:wfn1_param%nmo)
      double precision, intent(out) :: v_dcnb(1:wfn1_param%nmo)
c
      type(wfn1_error) :: error
c
      integer,          external :: wfn1_idmin
      type(wfn1_error), external :: wfn1_select_max
      type(wfn1_error), external :: wfn1_fetch_max
      type(wfn1_error), external :: wfn1_limit_max
c
      integer             npts        !< The number of points
      integer             mxpts       !< The maximum number of points
      parameter          (mxpts = 20) !< Set the number points
      double precision xx(mxpts)      !< The x-coordinates
      double precision fx(mxpts)      !< The function values
      double precision df(mxpts)      !< The function gradient
      double precision de             !< The energy gradient
      double precision tr             !< The trust radius
      double precision aa, bb, cc     !< The factor of the parabola
      double precision xm, fxm, dfm
      double precision tol
      double precision dnrm           !< Normalization factor
c
      integer indxmin
      integer ix1, ix2, ix3
      integer ipos          !< The position where xm was inserted
      logical o_eminimized  !< If b=0 then only redundant
                            !< transformations remain.
                            !< Interestingly even if a!=0 and b=0 the
                            !< only reasonable step size is 0. So b=0 is
                            !< the critical condition beyond which no
                            !< progress can be made minimizing the 
                            !< energy.
      logical o_reduced     !< Whether the gradient has been reduced
                            !< sufficiently
      logical o_limited     !< Have the Fock matrices been limited.
                            !< A limited Fock matrix has been reduced
                            !< to one element (the largest one from the
                            !< full Fock matrix).
c
      integer iproc         !< The current rank
      type(wfn1_en)  :: et  !< Temporary energy structure
      type(wfn1_wfn) :: wt  !< Temporary wave function structure
      type(wfn1_drv) :: dt  !< Temporary wave function derivative
                            !< structure
c
      character*17 pname
      parameter(pname = "wfn1_linesearch: ")
c
      iproc = ga_nodeid()
      tol = wfn1_param%convergence
      tr = 1.0d0
      o_limited = .false.
c
c     Store the initial point
c
 10   npts = 1
      xx       = 0.0d0
      fx(npts) = wfn1_energy%E_tot
      error = wfn1_etrace(wfn1_param,wfn1_wave,wfn1_deriv,wfn1_deriv)
      de = error%total
c     error = wfn1_trace(wfn1_deriv,wfn1_deriv)
      error = wfn1_select_max(wfn1_deriv)
      dnrm = error%total**2
      dnrm = 1.0d0/dnrm
      df(npts) = 1.0d0
c
c     Take one step
c
      npts = npts + 1
      xx(npts) = 1.0d0
      call wfn1_gen_trial_wfn(xx(npts),wfn1_param,wfn1_wave,wfn1_deriv,
     &                        wfn1_trial_wv,g_momo1,g_momo2)
      call wfn1_ortho(wfn1_mat,wfn1_trial_wv)
c     call wfn1_tst_ortho(wfn1_mat,wfn1_trial_wv,g_bfmo,g_momo1)
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &              wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &              g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &              g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &              v_dna,v_dnb,v_dca,v_dcb,
     &              v_dnna,v_dnnb,v_dcna,v_dcnb)
c     call wfn1_skew(wfn1_param,wfn1_trial_drv,
c    &               v_dna,v_dnb,v_dca,v_dcb)
c     call wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_trial_drv)
      fx(npts) = wfn1_trial_en%E_tot
c     error = wfn1_trace(wfn1_deriv,wfn1_trial_drv)
c     error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
      error = wfn1_fetch_max(wfn1_trial_drv)
      df(npts) = dnrm*(error%total**2)
      call wfn1_f0df0f1(xx(2),fx(1),de,fx(2),tr,aa,bb,cc,xm,fxm)
      if (xm.lt.-tr) then
        xm = -tr
      else if (xm.gt.1.0d0+tr) then
        xm =  1.0d0+tr
      endif
      if (abs(xm).le.tol) then
        if (xm.lt.0.0d0) then
          xm = -tol
        else
          xm =  tol
        endif
      endif
      call wfn1_fix_xm(npts,xx,xm)
c
c     Always select the next step (hoping that using the quadratic
c     approximation will pay off)
c
      call wfn1_gen_trial_wfn(xm,wfn1_param,wfn1_wave,wfn1_deriv,
     &                        wfn1_trial_wv,g_momo1,g_momo2)
      call wfn1_ortho(wfn1_mat,wfn1_trial_wv)
c     call wfn1_tst_ortho(wfn1_mat,wfn1_trial_wv,g_bfmo,g_momo1)
      call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &              wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &              g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &              g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &              v_dna,v_dnb,v_dca,v_dcb,
     &              v_dnna,v_dnnb,v_dcna,v_dcnb)
c     call wfn1_skew(wfn1_param,wfn1_trial_drv,
c    &               v_dna,v_dnb,v_dca,v_dcb)
c     error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
      error = wfn1_fetch_max(wfn1_trial_drv)
      dfm = dnrm*(error%total**2)
      call wfn1_insert(mxpts,npts,xx,fx,df,xm,wfn1_trial_en%E_tot,dfm,
     &                 ipos)
c
      o_eminimized = .false.
      o_reduced    = .false.
c
c     Minimize the energy
c
      do while ((.not.o_eminimized).and.(.not.o_reduced).and.
     &          npts.lt.mxpts)
c    & .or. 
c    &          wfn1_trial_en%E_tot > wfn1_energy%E_tot*(1.0d0+1.0d-8))
c
c       Predict a new xm
c
        indxmin = wfn1_idmin(npts,fx,1)
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
        call wfn1_f1f2f3(xx(1)-tr,xx(npts)+tr,
     &                   xx(ix1),xx(ix2),xx(ix3),
     &                   fx(ix1),fx(ix2),fx(ix3),
     &                   aa,bb,cc,xm,fxm)
        call wfn1_fix_xm(npts,xx,xm)
c
c       There are three different scenarios here:
c       1. INDXMIN.eq.1: Then XM must be XM < XX(2). If not then the
c          fitted parabola is useless and we set XM to the edge of the 
c          trust region XM = XX(1) - TR.
c       2. INDXMIN.eq.NPTS: Then XM must be XM > XX(NPTS-1). If not
c          then the fitted parabola is useless and we set XM to the 
c          edge of the trust region XM = XX(NPTS) + TR.
c       3. INDXMIN is anything else: Then XM must be XX(INDXMIN-1) < 
c          XM < XX(INDXMIN+1). If not then the fitted parabola is 
c          useless, instead we choose XM through bisecting the larger
c          of the intervals XX(INDXMIN-1):XX(INDXMIN) and 
c          XX(INDXMIN):XX(INDXMIN+1).
c
        if      (indxmin.eq.1.and.xm.ge.xx(2)) then
          xm = xx(1) - tr
        else if (indxmin.eq.npts.and.xm.le.xx(npts-1)) then
          xm = xx(npts) + tr
        else if (indxmin.gt.1.and.indxmin.lt.npts.and..not.
     &           (xx(indxmin-1).lt.xm.and.xm.lt.xx(indxmin+1))) then
          if (xx(indxmin)-xx(indxmin-1).gt.xx(indxmin+1)-xx(indxmin))
     &    then
            xm = 0.5d0*(xx(indxmin)+xx(indxmin-1))
          else
            xm = 0.5d0*(xx(indxmin+1)+xx(indxmin))
          endif
        endif
        if (abs(xm).le.tol) then
          if (xm.lt.0.0d0) then
            xm = -tol
          else
            xm =  tol
          endif
        endif
c
        call wfn1_gen_trial_wfn(xm,wfn1_param,wfn1_wave,wfn1_deriv,
     &                          wfn1_trial_wv,g_momo1,g_momo2)
        call wfn1_ortho(wfn1_mat,wfn1_trial_wv)
c       call wfn1_tst_ortho(wfn1_mat,wfn1_trial_wv,g_bfmo,g_momo1)
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                v_dna,v_dnb,v_dca,v_dcb,
     &                v_dnna,v_dnnb,v_dcna,v_dcnb)
c       call wfn1_skew(wfn1_param,wfn1_trial_drv,
c    &                 v_dna,v_dnb,v_dca,v_dcb)
c       call wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_trial_drv)
c       error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
        error = wfn1_fetch_max(wfn1_trial_drv)
        dfm = dnrm*(error%total**2)
        call wfn1_insert(mxpts,npts,xx,fx,df,xm,wfn1_trial_en%E_tot,dfm,
     &                   ipos)
c
c       Xm=xx(ipos) was predicted to be the minimum, test this
c       prediction
c
        call wfn1_closest_three(npts,xx,ipos,ix1,ix2,ix3)
        call wfn1_f1f2f3(xx(1)-tr,xx(npts)+tr,
     &                   xx(ix1),xx(ix2),xx(ix3),
     &                   fx(ix1),fx(ix2),fx(ix3),
     &                   aa,bb,cc,xm,fxm)
        if (abs(xm).le.tol) then
          if (xm.lt.0.0d0) then
            xm = -tol*100d0
          else
            xm =  tol*100d0
          endif
        endif
        o_reduced = (abs(2.0d0*aa*xx(ipos)+bb).lt.0.1d0*bb)
        o_eminimized = (abs(bb).lt.tol*0.01d0).and..not.o_reduced
      enddo
cHvD
      write(*,*)'HVD: ',npts,xm
      return
      if (wfn1_param%corr_expr.eq."wfn1s") return
cHvD
c     if (npts.ge.mxpts) then
c       write(*,*)"ipos   = ",ipos
c       write(*,*)"ix     = ",ix1,ix2,ix3
c       write(*,*)"a,b,c  = ",aa,bb,cc
c       write(*,*)"x,f(x) = ",xm,fxm
c       do ix1 = 1, mxpts
c         write(*,*)ix1,xx(ix1),fx(ix1)
c       enddo
c       indxmin = wfn1_idmin(npts,fx,1)
c       call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
c       call wfn1_f1f2f3(xx(1)-tr,xx(npts)+tr,
c    &                   xx(ix1),xx(ix2),xx(ix3),
c    &                   fx(ix1),fx(ix2),fx(ix3),
c    &                   aa,bb,cc,xm,fxm)
c       write(*,*)"indxmin = ",indxmin
c       write(*,*)"ix      = ",ix1,ix2,ix3
c       write(*,*)"a,b,c   = ",aa,bb,cc
c       write(*,*)"x,f(x)  = ",xm,fxm
c       call errquit("too many steps in linesearch",0,UERR)
c     endif
c
c     If the energy is minimized already and only redundant
c     transformations remain then minimize the norm of the Fock 
c     matrices.
c
c     This works only when minimizing one matrix element. So we need
c     to change the search direction to make this possible.
c     The change of search direction is effected by wfn1_limit_max.
c
      if (o_eminimized .and. npts.lt.mxpts .and. .not.o_limited) then
c
c       Limit the Fock matrix and redo the energy minimization on it
c       again.
c       The original Fock matrix provides a different search direction
c       than the limited Fock matrix. Hence we need to minimize the 
c       energy on the limited space first before attempting to minimize
c       the norm of the Fock matrix (under the assumption that the non-
c       zero elements correspond to redundant variables).
c
        wt = wfn1_wave
        wfn1_wave = wfn1_trial_wv
        wfn1_trial_wv = wt
c
        dt = wfn1_deriv
        wfn1_deriv = wfn1_trial_drv
        wfn1_trial_drv = dt
c
        wfn1_energy = wfn1_trial_en
c
        error = wfn1_select_max(wfn1_deriv)
        error = wfn1_limit_max(wfn1_deriv)
c
        o_reduced = .false.
        o_limited = .true.
c
        goto 10
      endif
      o_reduced = .false.
c     error = wfn1_limit_max(wfn1_deriv)
c
!     if (wfn1_param%corr_expr.eq."wfn1s") return
c
      do while (o_eminimized .and. npts.lt.mxpts .and. .not.o_reduced)
        indxmin = wfn1_idmin(npts,df,1)
        call wfn1_closest_three(npts,xx,indxmin,ix1,ix2,ix3)
        call wfn1_df1df2df3(xx(ix1),xx(ix2),xx(ix3),
     &                      df(ix1),df(ix2),df(ix3),
     &                      tr,aa,bb,cc,xm,fxm)
c       if (aa.lt.0.0d0) then
c         if (indxmin.eq.ix1) then
c           xm = 0.5d0*(xx(ix1)+xx(ix2))
c         else if (indxmin.eq.ix3) then
c           xm = 0.5d0*(xx(ix2)+xx(ix3))
c         else
c           call errquit(pname//" messed up",0,UERR)
c         endif
c       else
        call wfn1_fix_xm(npts,xx,xm)
c       endif
c
c       Shrink the trust region to avoid getting trapped generating
c       the same unsuitable answer over and over
c
c       tr = tr*0.95d0
c
        call wfn1_gen_trial_wfn(xm,wfn1_param,wfn1_wave,wfn1_deriv,
     &                          wfn1_trial_wv,g_momo1,g_momo2)
        call wfn1_ortho(wfn1_mat,wfn1_trial_wv)
c       call wfn1_tst_ortho(wfn1_mat,wfn1_trial_wv,g_bfmo,g_momo1)
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                v_dna,v_dnb,v_dca,v_dcb,
     &                v_dnna,v_dnnb,v_dcna,v_dcnb)
        call wfn1_skew(wfn1_param,wfn1_trial_drv,
     &                 v_dna,v_dnb,v_dca,v_dcb)
c       call wfn1_zero_vv(wfn1_param,wfn1_wave,wfn1_trial_drv)
c       error = wfn1_trace(wfn1_trial_drv,wfn1_trial_drv)
        error = wfn1_fetch_max(wfn1_trial_drv)
        dfm = dnrm*(error%total**2)
        call wfn1_insert(mxpts,npts,xx,fx,df,xm,wfn1_trial_en%E_tot,dfm,
     &                   ipos)
c
c       Test prediction
c
        call wfn1_closest_three(npts,xx,ipos,ix1,ix2,ix3)
        call wfn1_df1df2df3(xx(ix1),xx(ix2),xx(ix3),
     &                      df(ix1),df(ix2),df(ix3),
     &                      tr,aa,bb,cc,xm,fxm)
        o_reduced = (abs((2.0d0*aa*xx(ipos)+bb)/bb).lt.0.1d0)
      enddo
c
      end subroutine wfn1_linesearch
c
c-----------------------------------------------------------------------
