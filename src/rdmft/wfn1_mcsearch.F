c-----------------------------------------------------------------------
c
C>
C> \brief Perform a Monte-Carlo search towards the minimum total energy
C>
C> Given a starting wave function and energy we generate an orthogonal
C> matrix, and apply this matrix as a 1-particle function update.
C> Dependent on the effect of the update on the energy we accept the
C> updated wave function with a particular probability.
C>
C> The orthogonal transformation is generated from a unit matrix,
C> where the off-diagonal elements are replaced with random numbers 
C> in the range [-a,a]. The values "a" will have to be relatively small
C> wrt. 1.0. This is required so that we can preserve the order of the
C> 1-particle functions in the updated wave function. After we have
C> filled this matrix we use modified Gramm-Schmidt to orthonormalize
C> the columns.
C>
C> The Monte Carlo approach involves guessing a step. First the
C> forward step is tried, if that does not lower the energy then
C> the backward step is tried, if that also does not lower the 
C> energy one has enough information to try a step that corresponds to
C> the minimum of a quadratic approximation. If that also fails then
C> simply try a new random step.
C>
C> Assume that the energy is a quadratic function
C> \f{eqnarray*}{
C>   E &=& ax^2 + bx + c
C> \f}
C> then with the three points calculated we have
C> \f{eqnarray*}{
C>   E_{-1} &=& a - b + c \\
C>   E_{ 0} &=& c \\
C>   E_{ 1} &=& a + b + c \\
C>   E_{1} + E_{-1} &=& 2a + 2c \\
C>   E_{1} - E_{-1} &=& 2b \\
C>   c &=& E_{0} \\
C>   b &=& (E_{1} - E_{-1})/2 \\
C>   a &=& (E_{1} + E_{-1} - 2E_{0})/2
C> \f}
C> Minimizing the energy gives
C> \f{eqnarray*}{
C>   \frac{\partial E}{\partial x} &=& 2ax + b = 0 \\
C>   x &=& -\frac{b}{2a}
C> \f}
C>
      subroutine wfn1_mcsearch(wfn1_param,wfn1_inst,wfn1_mat,
     &                         wfn1_wave,
!    &                         wfn1_deriv,
     &                         wfn1_energy,wfn1_trial_wv,
!    &                         wfn1_trial_drv,wfn1_part_drv,
     &                         wfn1_trial_en,
     &                         wfn1_mcbox,wfn1_rejct,wfn1_accpt,
     &                         g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                         g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                         g_momo5,g_momo6,
     &                         g_4t1,g_4t2,
     &                         gz_momo1,gz_momo2,gz_momo3,gz_momo4,
     &                         fock_count)
!    &                         v_dna,v_dnb,v_dca,v_dcb,
!    &                         v_dnna,v_dnnb,v_dcna,v_dcnb)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
!include "wfn1_trace.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
!include "wfn1_en.fh"
c
      type(wfn1_prm),     intent(inout) :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(inout) :: wfn1_wave
!     type(wfn1_drv),     intent(inout) :: wfn1_deriv
      type(wfn1_en),      intent(inout) :: wfn1_energy
      type(wfn1_wfn),     intent(inout) :: wfn1_trial_wv
!     type(wfn1_drv),     intent(inout) :: wfn1_trial_drv
!     type(wfn1_drv),     intent(inout) :: wfn1_part_drv
      type(wfn1_en),      intent(inout) :: wfn1_trial_en
      type(wfn1_bx),      intent(inout) :: wfn1_mcbox
      type(wfn1_bx),      intent(inout) :: wfn1_rejct
      type(wfn1_bx),      intent(inout) :: wfn1_accpt
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
      integer, intent(inout) :: g_momo5 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_momo6 ! nmo x nmo scratch GA
      integer, intent(inout) :: g_4t1   ! nbf x nbf x nbf x nbf scratch GA
      integer, intent(inout) :: g_4t2   ! nbf x nbf x nbf x nbf scratch GA
      integer, intent(inout) :: gz_momo1 ! nmo x nmo complex scratch GA
      integer, intent(inout) :: gz_momo2 ! nmo x nmo complex scratch GA
      integer, intent(inout) :: gz_momo3 ! nmo x nmo complex scratch GA
      integer, intent(inout) :: gz_momo4 ! nmo x nmo complex scratch GA
c
      integer, intent(inout) :: fock_count ! # energy evaluations
c
c     Natural orbital energies
c
!     double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
c
c     Correlation function energies
c
!     double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
c
c     Natural orbital expectation values of the Correlation
c     function Fock matrices
c
!     double precision, intent(out) :: v_dnna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dnnb(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcna(1:wfn1_param%nmo)
!     double precision, intent(out) :: v_dcnb(1:wfn1_param%nmo)
c
      integer :: iproc ! the current rank
      integer :: nproc ! the number of ranks
      integer :: ilo, ihi, jlo, jhi ! block limits
      integer :: ii ! counter
      integer :: idum ! dummy
      integer :: nmo  ! number of MOs
      integer :: ntries ! number of attempts
c
      double precision, allocatable :: buf(:,:)
      double precision :: rlim ! random number limits
      logical :: oclosed  ! is this a closed shell calculation?
      logical :: ostop    ! do we stop iterating?
      logical :: okeep    ! do we keep the quadratic approximation
      double precision :: Ewf       ! Current wave function energy
      double precision :: Ewf_trial ! Trial wave function energy
      double precision :: random    ! One random number
      double precision :: de        ! Change in energy
      double precision :: pfac      ! Probability scale factor
      double precision :: sfac      ! Dimensionality adjustment factor
      double precision :: tfac      ! Dimensionality adjustment factor
      double precision :: Em1       ! \f$E_{-1}\f$
      double precision :: E0        ! \f$E_{ 0}\f$
      double precision :: Emin      ! \f$E_{min}\f$ predicted minimum
      double precision :: Ep1       ! \f$E_{+1}\f$
      double precision :: aa
      double precision :: bb
      double precision :: cc
      double precision :: xx
c
      iproc   = ga_nodeid()
      nproc   = ga_nnodes()
      oclosed = wfn1_param%nea.eq.wfn1_param%neb
      nmo     = wfn1_param%nmo
      rlim    = wfn1_param%rlim
      pfac    = wfn1_param%pfac
      Ewf     = wfn1_energy%E_tot
      E0      = Ewf
c     sfac    = 4*(nmo*(nmo-1)/2)
c     sfac    = 1.0d0/sfac 
c     tfac    = 1.08d0**sfac
c     sfac    = 0.5d0**sfac
!
!     sfac = shrink factor which we choose independent of the number
!            of dimensions in the problem
!     tfac = was the growth factor that was used to accelerate the 
!            calculation upon successful steps
!
      sfac = 0.5d0
      tfac = 1.08d0
!
      ntries  = 0
!     call ga_distribution(g_momo1,iproc,ilo,ihi,jlo,jhi)
!     allocate(buf(ilo:ihi,jlo:jhi))
      call wfn1_init_box(wfn1_rejct,0.0d0)
      call wfn1_init_box(wfn1_accpt,0.0d0)
      do
        ntries = ntries + 1
cDEBUG
c       write(*,*)ntries,fock_count,E0,rlim
cDEBUG
c
c       Natural orbitals random
c
        call wfn1_gen_random_skew(g_momo1,wfn1_mcbox%g_na,rlim)
        if (oclosed) then
          call ga_copy(g_momo1,g_momo2)
        else
          call wfn1_gen_random_skew(g_momo2,wfn1_mcbox%g_nb,rlim)
        endif
c
c       Correlation functions random
c
        call wfn1_gen_random_skew(gz_momo1,wfn1_mcbox%g_ca,rlim)
        if (oclosed) then
          call ga_copy(gz_momo1,gz_momo2)
        else
          call wfn1_gen_random_skew(gz_momo2,wfn1_mcbox%g_cb,rlim)
        endif
c
c       Select an direction orthogonal to previously minimized ones
c
        call wfn1_dir_orthog_vec(g_momo1,g_momo2,gz_momo1,gz_momo2)
c
c       Create rotations
c
        call wfn1_skew_to_rot(g_momo1,g_momo3)
        call wfn1_skew_to_rot(g_momo2,g_momo4)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,
     &                 g_momo3,0.0d0,wfn1_trial_wv%g_na)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,
     &                 g_momo4,0.0d0,wfn1_trial_wv%g_nb)
c
        call wfn1_skew_to_rot(gz_momo1,gz_momo3)
        call wfn1_skew_to_rot(gz_momo2,gz_momo4)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &                 gz_momo3,0.0d0,wfn1_trial_wv%g_ca)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &                 gz_momo4,0.0d0,wfn1_trial_wv%g_cb)
c
c       Evaluate energy
c
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
!    &                wfn1_trial_drv,wfn1_part_drv,
     &                wfn1_trial_en,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,g_momo3,g_momo4,g_momo5,g_momo6,
     &                g_4t1,g_4t2)
!    &                v_dna,v_dnb,v_dca,v_dcb,
!    &                v_dnna,v_dnnb,v_dcna,v_dcnb)
        fock_count = fock_count + 1

        Ewf_trial = wfn1_trial_en%E_tot
        Em1       = Ewf_trial
        de = Ewf_trial - Ewf
        if (iproc.eq.0) then
c
c         If the energy goes down we are going in the right direction
c
          ostop = de.lt.0.0d0
        endif
        call ga_brdcst(1200,ostop,ma_sizeof(MT_LOG,1,MT_BYTE),0)
        if (.not.ostop) then
cDEBUG
c         write(*,*)"reject 1"
c         call util_flush(6)
cDEBUG
          call wfn1_sgn_offdiag(g_momo1)
          call wfn1_sgn_offdiag(g_momo2)
          call wfn1_skew_to_rot(g_momo1,g_momo3)
          call wfn1_skew_to_rot(g_momo2,g_momo4)
          call wfn1_sgn_offdiag(gz_momo1)
          call wfn1_sgn_offdiag(gz_momo2)
          call wfn1_skew_to_rot(gz_momo1,gz_momo3)
          call wfn1_skew_to_rot(gz_momo2,gz_momo4)
          call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,
     &                   g_momo3,0.0d0,wfn1_trial_wv%g_na)
          call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,
     &                   g_momo4,0.0d0,wfn1_trial_wv%g_nb)
          call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &                   gz_momo3,0.0d0,wfn1_trial_wv%g_ca)
          call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &                   gz_momo4,0.0d0,wfn1_trial_wv%g_cb)
          call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                  wfn1_trial_en,
     &                  g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                  g_bfmo,g_momo3,g_momo4,g_momo5,g_momo6,
     &                  g_4t1,g_4t2)
          fock_count = fock_count + 1
          Ewf_trial  = wfn1_trial_en%E_tot
          Ep1        = Ewf_trial
          de = Ewf_trial - Ewf
          if (iproc.eq.0) then
            ostop = de.lt.0.0d0
          endif
          call ga_brdcst(1300,ostop,ma_sizeof(MT_LOG,1,MT_BYTE),0)
        endif
        if (.not.ostop) then
cDEBUG
c         write(*,*)"reject 2", Em1,E0,Ep1
c         call util_flush(6)
cDEBUG
c         The raw step was rejected so update the reject box
          call wfn1_max_box(wfn1_rejct,g_momo1,g_momo2,
     &                      gz_momo1,gz_momo2)
          cc   = E0
          bb   = (Ep1-Em1)/2.0d0
          aa   = (Ep1+Em1-2.0d0*E0)/2.0d0
          if (abs(aa).gt.1.0d-10) then
            xx = -bb/(2.0d0*aa)
          else
            xx = 0.0d0
          endif
          Emin = aa*xx*xx+bb*xx+cc
          if (abs(xx).gt.1.0d-5) then
c           call wfn1_write_parabola("parabole",fock_count,aa,bb,cc)
c           call wfn1_write_line("line",fock_count+1,
c    &           wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
c    &           wfn1_trial_wv,wfn1_trial_en,
c    &           g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
c    &           g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,g_momo5,g_momo6,
c    &           gz_momo1,gz_momo2,gz_momo3,gz_momo4,
c    &           g_4t1,g_4t2)
            call ga_scale(g_momo1,xx)
            call ga_scale(g_momo2,xx)
            call wfn1_skew_to_rot(g_momo1,g_momo3)
            call wfn1_skew_to_rot(g_momo2,g_momo4)
            call ga_scale(gz_momo1,xx)
            call ga_scale(gz_momo2,xx)
            call wfn1_skew_to_rot(gz_momo1,gz_momo3)
            call wfn1_skew_to_rot(gz_momo2,gz_momo4)
            call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,
     &                     g_momo3,0.0d0,wfn1_trial_wv%g_na)
            call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,
     &                     g_momo4,0.0d0,wfn1_trial_wv%g_nb)
            call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &                     gz_momo3,0.0d0,wfn1_trial_wv%g_ca)
            call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &                     gz_momo4,0.0d0,wfn1_trial_wv%g_cb)
            call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                    wfn1_trial_en,
     &                    g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                    g_bfmo,g_momo3,g_momo4,g_momo5,g_momo6,
     &                    g_4t1,g_4t2)
            fock_count = fock_count + 1
            Ewf_trial  = wfn1_trial_en%E_tot
            de = Ewf_trial - Ewf
            if (iproc.eq.0) then
              okeep = de.lt.0.0d0
cDEBUG
c             write(*,*)'HVD 0',fock_count,Em1,E0,Ep1
c             write(*,*)'HVD 1',fock_count,Emin,Ewf_trial,Ewf_trial-Emin
cDEBUG
            endif
            call ga_brdcst(1300,okeep,ma_sizeof(MT_LOG,1,MT_BYTE),0)
            if (okeep) then
              call wfn1_copy_wfn(wfn1_trial_wv,wfn1_wave)
              Ewf = wfn1_trial_en%E_tot
              call wfn1_max_box(wfn1_accpt,g_momo1,g_momo2,
     &                          gz_momo1,gz_momo2)
              call wfn1_dir_store_vec()
c             call wfn1_write_random("direction",fock_count,g_momo1,
c    &                               g_momo2,gz_momo1,gz_momo2)
c             call wfn1_write_points("points",fock_count,
c    &                               -1.0d0,0.0d0,1.0d0,xx,
c    &                               Em1,E0,Ep1,Emin,Ewf)
              E0  = Ewf
            endif
cDEBUG
c           if (okeep) then
c             write(*,*)'accept 3',xx
c           else
c             write(*,*)'reject 3',xx
c           endif
c           call util_flush(6)
cDEBUG
          endif
        endif
        if (.not.ostop) then
!
!         The step was rejected so update the rejects box
!
!         This is now recorded before the quadratic update step.
!
c         call wfn1_max_box(wfn1_rejct,g_momo1,g_momo2,
c    &                      gz_momo1,gz_momo2)
        else
!
!         The step was accepted so update the accepts box
!
          call wfn1_max_box(wfn1_accpt,g_momo1,g_momo2,
     &                      gz_momo1,gz_momo2)
!
!         Reshape the Monte Carlo box
!
          call wfn1_reshape_box(wfn1_param,wfn1_mcbox,wfn1_rejct,
     &                          wfn1_accpt)
        endif
        if (.not.ostop.and.(mod(ntries,20).eq.0)) then
c         call wfn1_resquare_box(wfn1_param,wfn1_mcbox)
!
!         Check whether we have reached the convergence criteria
!         already. If so, quit the calculation (no matter what).
!
          ostop = 2.0d0*wfn1_param%rlim*tfac.lt.wfn1_param%convergence
          if (ostop) then
            call wfn1_copy_wfn(wfn1_wave,wfn1_trial_wv)
          endif
          wfn1_param%rlim = sfac*wfn1_param%rlim
          rlim = wfn1_param%rlim
        endif
        if (ostop) then
c         accept the trial wave function
          wfn1_param%rlim = tfac*wfn1_param%rlim
          exit
        endif
      end do
!     deallocate(buf)
c
      end
c
c-----------------------------------------------------------------------
c
C> \brief Generate a random skew symmetric matrix
C>
C> We will use this matrix to set a rotation matrix up. For every 2x2
C> rotation we have that
C>
C>      (R_ii R_ij)   ( c  -s )
C>  R = (         ) = (       )
C>      (R_ji R_jj)   ( s   c )
C>
C>  where R_ij = -R_ji (see e.g.
C>  https://en.wikipedia.org/wiki/Jacobi_rotation). 
C>  If we have a skew symmetric matrix T then we have
C>  * R_ii = R_jj = 1/sqrt(T_ij**2+1)
C>  * R_ij = T_ij/sqrt(T_ij**2+1)
C>  * R_ji = T_ji/sqrt(T_ji**2+1)
C>  For many 2x2 rotations R'_ii = \prod_j 1/sqrt(T_ij**2+1)
C>
      subroutine wfn1_gen_random_skew(g_rand,g_box,rlim)
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_rand ! Global array
      integer, intent(in) :: g_box  ! Global array
      double precision, intent(in)  :: rlim ! random number limits
c
      double precision, allocatable :: buf(:,:,:)
      double precision, allocatable :: buft(:,:,:)
      complex(kind=8),  allocatable :: zbuf(:,:)
      complex(kind=8),  allocatable :: zbuft(:,:)
      double precision, allocatable :: box(:,:)
c
      integer :: idum
      integer :: ii ! counter
      integer :: ilo, ihi, jlo, jhi
      integer :: iproc
      integer :: itype, jj
      double precision :: tt
c
      iproc = ga_nodeid()
      call ga_inquire(g_rand,itype,ii,jj)
      call ga_distribution(g_rand,iproc,ilo,ihi,jlo,jhi)
      call ga_zero(g_rand)
      if (itype.eq.MT_DBL) then
        if (ilo.ge.jlo) then
          allocate(buf(1:1,ilo:ihi,jlo:jhi))
          allocate(buft(1:1,jlo:jhi,ilo:ihi))
          allocate(box(ilo:ihi,jlo:jhi))
          call random_number(buf)
          call ga_get(g_box,ilo,ihi,jlo,jhi,box,ihi-ilo+1)
          do jj = jlo, jhi
            do ii = ilo, ihi
              if (ii.lt.jj) then
                buf(1,ii,jj)  = (2.0d0*buf(1,ii,jj)-1.0d0)
     &                        * rlim*box(ii,jj)
                buft(1,jj,ii) = -buf(1,ii,jj)
              else if (ii.gt.jj) then
                buf(1,ii,jj)  = 0.0d0
                buft(1,jj,ii) = 0.0d0
              else
                buf(1,ii,jj)  = 0.0d0
                buft(1,jj,ii) = 0.0d0
              endif
            enddo
          enddo
          call ga_acc(g_rand,ilo,ihi,jlo,jhi,buf, ihi-ilo+1,1.0d0)
          call ga_acc(g_rand,jlo,jhi,ilo,ihi,buft,jhi-jlo+1,1.0d0)
          deallocate(box)
          deallocate(buft)
          deallocate(buf)
        endif
      else if (itype.eq.MT_DCPL) then
        if (ilo.ge.jlo) then
          allocate(buf(1:2,ilo:ihi,jlo:jhi))
          allocate(zbuf(ilo:ihi,jlo:jhi))
          allocate(zbuft(jlo:jhi,ilo:ihi))
          allocate(box(ilo:ihi,jlo:jhi))
          call random_number(buf)
          call ga_get(g_box,ilo,ihi,jlo,jhi,box,ihi-ilo+1)
          buf = (2.0d0*buf-1.0d0)*rlim
          do jj = jlo, jhi
            do ii = ilo, ihi
              if (ii.lt.jj) then
                zbuf(ii,jj)  = complex(buf(1,ii,jj),buf(2,ii,jj))
     &                       * box(ii,jj)
                zbuft(jj,ii) = -zbuf(ii,jj)
              else if (ii.gt.jj) then
                zbuf(ii,jj)  = complex(0.0d0,0.0d0)
                zbuft(jj,ii) = complex(0.0d0,0.0d0)
              else
                zbuf(ii,jj)  = complex(0.0d0,0.0d0)
                zbuft(jj,ii) = complex(0.0d0,0.0d0)
              endif
            enddo
          enddo
          call ga_acc(g_rand,ilo,ihi,jlo,jhi,zbuf, ihi-ilo+1,
     &                complex(1.0d0,0.0d0))
          call ga_acc(g_rand,jlo,jhi,ilo,ihi,zbuft,jhi-jlo+1,
     &                complex(1.0d0,0.0d0))
          deallocate(box)
          deallocate(zbuft)
          deallocate(zbuf)
          deallocate(buf)
        endif
      endif
      call ga_sync()
c
      end
c
c-----------------------------------------------------------------------
c
!> \brief form a rotation matrix from a skew symmetric matrix
!>
      subroutine wfn1_skew_to_rot(g_skew,g_rot)
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_skew ! Global array
      integer, intent(in) :: g_rot  ! Global array
      integer             :: ioption ! choose the transformation:
                                     ! 1. tangens based
                                     ! 2. angle based
c
      double precision, allocatable :: buf(:,:)
      complex(kind=8),  allocatable :: zbuf(:,:)
c
      integer :: iidim, jjdim
      integer :: ii ! counter
      integer :: ilo, ihi, jlo, jhi
      integer :: ncol
      integer :: iproc
      integer :: nproc
      integer :: itype, jj
      double precision :: cc
      double precision :: ct
c
      ioption = 2
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      call ga_inquire(g_skew,itype,iidim,jjdim)
      ilo  = 1
      ihi  = iidim
      ncol = (jjdim+nproc-1)/nproc
      jlo  = iproc*ncol+1
      jhi  = min((iproc+1)*ncol,jjdim)
      if (itype.eq.MT_DBL) then
        if (jlo.le.jhi) then
          allocate(buf(ilo:ihi,jlo:jhi))
          call ga_get(g_skew,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
          do jj = jlo, jhi
            ct = 1.0d0
            do ii = ilo, ihi
              if (ioption.eq.1) then
                cc = 1.0d0/sqrt(buf(ii,jj)*buf(ii,jj)+1.0d0)
              elseif (ioption.eq.2) then
                cc = cos(buf(ii,jj))
              endif
              ct = ct*cc
              if (ioption.eq.1) then
                buf(ii,jj) = buf(ii,jj)*cc
              elseif (ioption.eq.2) then
                buf(ii,jj) = sin(buf(ii,jj))
              endif
            enddo
            buf(jj,jj) = ct
          enddo
          call ga_put(g_rot,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
          deallocate(buf)
        endif
      else if (itype.eq.MT_DCPL) then
        if (jlo.le.jhi) then
          allocate(zbuf(ilo:ihi,jlo:jhi))
          call ga_get(g_skew,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
          do jj = jlo, jhi
            ct = 1.0d0
            do ii = ilo, ihi
              if (ioption.eq.1) then
                cc = 1.0d0/sqrt(abs(zbuf(ii,jj))**2+1.0d0)
              elseif (ioption.eq.2) then
                cc = cos(abs(zbuf(ii,jj)))
              endif
              ct = ct*cc
              if (ioption.eq.1) then
                zbuf(ii,jj) = zbuf(ii,jj)*cc
              elseif (ioption.eq.2) then
                zbuf(ii,jj) = sin(zbuf(ii,jj))
              endif
            enddo
            zbuf(jj,jj) = complex(ct,0.0d0)
          enddo
          call ga_put(g_rot,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
          deallocate(zbuf)
        endif
      endif
      call ga_sync()
      call wfn1_orthog_vecs(g_rot,iidim,.false.)
c
      end
c
c-----------------------------------------------------------------------
