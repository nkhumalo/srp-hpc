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
      subroutine wfn1_mcsearch(wfn1_param,wfn1_inst,wfn1_mat,
     &                         wfn1_wave,wfn1_deriv,wfn1_energy,
     &                         wfn1_trial_wv,wfn1_trial_drv,
     &                         wfn1_part_drv,wfn1_trial_en,
     &                         g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                         g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                         v_dna,v_dnb,v_dca,v_dcb,
     &                         v_dnna,v_dnnb,v_dcna,v_dcnb)
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
      double precision :: Ewf       ! Current wave function energy
      double precision :: Ewf_trial ! Trial wave function energy
      double precision :: random    ! One random number
      double precision :: de        ! Change in energy
      double precision :: pfac      ! Probability scale factor
      double precision :: sfac      ! Dimensionality adjustment factor
c
      iproc   = ga_nodeid()
      nproc   = ga_nnodes()
      oclosed = wfn1_param%nea.eq.wfn1_param%neb
      nmo     = wfn1_param%nmo
      rlim    = wfn1_param%rlim
      pfac    = wfn1_param%pfac
      Ewf     = wfn1_energy%E_tot
      sfac    = wfn1_param%nea*(nmo-wfn1_param%nea)
      sfac    = 1.0d0/sfac 
      sfac    = 0.5d0**sfac
      ntries  = 0
!     call ga_distribution(g_momo1,iproc,ilo,ihi,jlo,jhi)
!     allocate(buf(ilo:ihi,jlo:jhi))
      do
        ntries = ntries + 1
c
c       Natural orbitals
c
        call wfn1_gen_random_orthog(g_momo1,rlim)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_na,
     &                 g_momo1,0.0d0,wfn1_trial_wv%g_na)
        if (.not.oclosed) then
          call wfn1_gen_random_orthog(g_momo1,rlim)
        endif
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_nb,
     &                 g_momo1,0.0d0,wfn1_trial_wv%g_nb)
c
c       Correlation functions
c
        call wfn1_gen_random_orthog(g_momo1,rlim)
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
     &                 g_momo1,0.0d0,wfn1_trial_wv%g_ca)
        if (.not.oclosed) then
          call wfn1_gen_random_orthog(g_momo1,rlim)
        endif
        call wfn1_gemm('n','n',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
     &                 g_momo1,0.0d0,wfn1_trial_wv%g_cb)
c
c       Evaluate energy
c
        call wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_trial_wv,
     &                wfn1_trial_drv,wfn1_part_drv,wfn1_trial_en,
     &                g_bfbf1,g_bfbf2,g_bfbf3,g_bfbf4,
     &                g_bfmo,g_momo1,g_momo2,g_momo3,g_momo4,
     &                v_dna,v_dnb,v_dca,v_dcb,
     &                v_dnna,v_dnnb,v_dcna,v_dcnb)

        Ewf_trial = wfn1_trial_en%E_tot
        de = Ewf - Ewf_trial
        if (iproc.eq.0) then
c         call random_number(random)
c         ostop = random.le.exp(pfac*de)
          ostop = de.gt.0.0d0
        endif
        call ga_brdcst(1200,ostop,ma_sizeof(MT_LOG,1,MT_BYTE),0)
        if (.not.ostop.and.(mod(ntries,100).eq.0)) then
c         ostop = .true.
          wfn1_param%rlim = sfac*wfn1_param%rlim
          rlim = wfn1_param%rlim
c         wfn1_trial_en = wfn1_energy
c         call ga_copy(wfn1_wave%g_na,wfn1_trial_wv%g_na)
c         call ga_copy(wfn1_wave%g_nb,wfn1_trial_wv%g_nb)
c         call ga_copy(wfn1_wave%g_ca,wfn1_trial_wv%g_ca)
c         call ga_copy(wfn1_wave%g_cb,wfn1_trial_wv%g_cb)
c         call ga_copy(wfn1_deriv%g_dna,wfn1_trial_drv%g_dna)
c         call ga_copy(wfn1_deriv%g_dnb,wfn1_trial_drv%g_dnb)
c         call ga_copy(wfn1_deriv%g_dca,wfn1_trial_drv%g_dca)
c         call ga_copy(wfn1_deriv%g_dcb,wfn1_trial_drv%g_dcb)
        endif
        if (ostop) then
c         accept the trial wave function
cDEBUG
          if (iproc.eq.0) then
          write(*,*)'HVD: mc: ',ntries,Ewf,Ewf_trial
          endif
cDEBUG
          exit
        endif
      end do
!     deallocate(buf)
c
      end
c
c-----------------------------------------------------------------------
c
C> \brief Generate a random orthogonal matrix
C>
      subroutine wfn1_gen_random_orthog(g_rand,rlim)
      implicit none
c
#include "global.fh"
#include "mafdecls.fh"
c
      integer, intent(in) :: g_rand ! Global array
      double precision, intent(in)  :: rlim ! random number limits
c
      double precision, allocatable :: buf(:,:,:)
      complex(kind=8),  allocatable :: zbuf(:,:)
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
      if (itype.eq.MT_DBL) then
        allocate(buf(1:1,ilo:ihi,jlo:jhi))
        call random_number(buf)
        buf = (2.0d0*buf-1.0d0)*rlim
        do ii = max(ilo,jlo),min(ihi,jhi)
          buf(1,ii,ii) = 1.0d0
        enddo
        call ga_put(g_rand,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
        deallocate(buf)
      else if (itype.eq.MT_DCPL) then
        allocate(buf(1:2,ilo:ihi,jlo:jhi))
        allocate(zbuf(ilo:ihi,jlo:jhi))
        call random_number(buf)
        buf = (2.0d0*buf-1.0d0)*rlim
        do jj = jlo, jhi
          do ii = ilo, ihi
            zbuf(ii,jj) = complex(buf(1,ii,jj),buf(2,ii,jj))
          enddo
        enddo
        do ii = max(ilo,jlo),min(ihi,jhi)
          zbuf(ii,ii) = complex(1.0d0,0.0d0)
        enddo
        call ga_put(g_rand,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
        deallocate(zbuf)
        deallocate(buf)
      endif
      call ga_sync()
      call wfn1_orthog(g_rand,idum,.false.)
c
      end
c
c-----------------------------------------------------------------------
