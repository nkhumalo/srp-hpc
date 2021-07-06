c-----------------------------------------------------------------------
c
C> \brief Calculate the energy and the derivative of the energy wrt
C> the wavefunction parameters
C>
C> This routine calculates the various energy contributions. At the
C> moment only a Monte Carlo optimizer is supported for which no
C> derivatives are needed.
C>
      subroutine wfn1_e_g(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                    wfn1_energy,
     &                    g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,
     &                    g_bfmo,
     &                    g_momo,g_momo2,g_momo3,g_momo4,
     &                    g_4t1,g_4t2)
      implicit none
C 
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
!include "wfn1_matcnst.fh"
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
!include "wfn1_en.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
!     type(wfn1_drv),     intent(inout) :: wfn1_deriv
!     type(wfn1_drv),     intent(inout) :: wfn1_part_drv
      type(wfn1_en),      intent(out)   :: wfn1_energy
c
      integer, intent(in) :: g_bfbf  ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf2 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf3 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf4 ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfmo  ! nbf x nmo scratch GA
      integer, intent(in) :: g_momo  ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo2 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo3 ! nmo x nmo scratch GA
      integer, intent(in) :: g_momo4 ! nmo x nmo scratch GA
      integer, intent(in) :: g_4t1   ! nbf x nbf x nbf x nbf scratch GA
      integer, intent(in) :: g_4t2   ! nbf x nbf x nbf x nbf scratch GA
c
!>    The \$\alpha\f$ natural orbital energies
!     double precision, intent(out) :: v_dna(1:wfn1_param%nmo)
!>    The \$\beta\f$ natural orbital energies
!     double precision, intent(out) :: v_dnb(1:wfn1_param%nmo)
!>    The \$\alpha\f$ correlation function energies
!     double precision, intent(out) :: v_dca(1:wfn1_param%nmo)
!>    The \$\beta\f$ correlation function energies
!     double precision, intent(out) :: v_dcb(1:wfn1_param%nmo)
!>    The \$\alpha\f$ natural orbital expectation values of
!>    natural orbital Fock matrix
!     double precision, intent(out) :: v_dnna(1:wfn1_param%nmo)
!>    The \$\beta\f$ natural orbital expectation values of
!>    natural orbital Fock matrix
!     double precision, intent(out) :: v_dnnb(1:wfn1_param%nmo)
!>    The \$\alpha\f$ natural orbital expectation values of
!>    correlation function Fock matrix
!     double precision, intent(out) :: v_dcna(1:wfn1_param%nmo)
!>    The \$\beta\f$ natural orbital expectation values of
!>    correlation function Fock matrix
!     double precision, intent(out) :: v_dcnb(1:wfn1_param%nmo)
!>    Non-correlation, non-thermal expectation values
!     double precision, allocatable :: v_dpa(:)
!     double precision, allocatable :: v_dpb(:)
!>
      double precision, allocatable :: buf(:,:)
      double precision, allocatable :: bufji(:,:)
      complex(kind=8),  allocatable :: zbuf(:,:)
c
      integer ilo,  ihi  ! range limits
      integer jlo,  jhi  ! range limits
      integer mxlo, mnhi ! range limits
      integer ii         ! counter
      integer jj         ! counter
      integer nmo        ! local copy of the number of MOs
c
      integer nproc
      integer iproc
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
c     Zero all relevant arrays so we can just add contributions
c
c     call ga_zero(g_bfmo)
c     call ga_zero(g_momo)
!     call ga_zero(wfn1_deriv%g_dna)
!     call ga_zero(wfn1_deriv%g_dnb)
!     call ga_zero(wfn1_deriv%g_dca)
!     call ga_zero(wfn1_deriv%g_dcb)
!     call ga_zero(wfn1_part_drv%g_dna)
!     call ga_zero(wfn1_part_drv%g_dnb)
!     call ga_zero(wfn1_part_drv%g_dca)
!     call ga_zero(wfn1_part_drv%g_dcb)
      call ga_sync()
cDEBUG
c     call wfn1_print_max(wfn1_deriv,"wfn1_e_g A")
c     write(*,*)'HVD: wfn1_e_g'
c     call wfn1_print_wave(wfn1_wave)
cDEBUG
c
      wfn1_energy%E_tot    = 0.0d0
      wfn1_energy%E_1el    = 0.0d0
      wfn1_energy%E_2el    = 0.0d0
      wfn1_energy%E_corr   = 0.0d0
      wfn1_energy%E_nuc    = wfn1_param%Enuc
      wfn1_energy%E_ent    = 0.0d0
      wfn1_energy%E_dft    = 0.0d0
c
      nmo = wfn1_param%nmo
c
c     Arrange memory for the diagonal correlation function Fock
c     matrices in natural orbital basis
c
!     allocate(v_dpa(1:wfn1_param%nmo))
!     allocate(v_dpb(1:wfn1_param%nmo))
c
!     v_dna  = 0.0d0
!     v_dnb  = 0.0d0
!     v_dca  = 0.0d0
!     v_dcb  = 0.0d0
!     v_dpa  = 0.0d0
!     v_dpb  = 0.0d0
!     v_dnna = 0.0d0
!     v_dnnb = 0.0d0
!     v_dcna = 0.0d0
!     v_dcnb = 0.0d0
c
c     Evaluate the various Fock matrix terms
c     --------------------------------------
c
c     The 2-electron terms
c
      if (wfn1_param%corr_expr == "wfn1_mx") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
      else if (wfn1_param%corr_expr == "wfn1s") then
        call wfn1_e_g_2el_s(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                      wfn1_deriv,wfn1_part_drv,v_dca,v_dcb,
!    &                      v_dpa,v_dpb,
     &                      g_momo,g_momo2,
!    &                      g_momo3,g_momo4,
     &                      wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
      else if (wfn1_param%corr_expr == "wfn1c1") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                         wfn1_part_drv,v_dpa,v_dpb,
     &                         g_momo,g_momo2,
!    &                         g_momo3,g_momo4,
     &                         wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c2") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr2(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c3") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr3(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c4") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr4(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c5") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr5(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c6") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr6(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c7") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr7(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c8") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr8(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c9") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                       wfn1_part_drv,v_dpa,v_dpb,
     &                       g_momo,g_momo2,
!    &                       g_momo3,g_momo4,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr9(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                          wfn1_part_drv,v_dpa,v_dpb,
     &                          g_momo,g_momo2,
!    &                          g_momo3,g_momo4,
     &                          wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c10") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       g_momo,g_momo2,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr10(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,
     &                           g_momo,g_momo2,
     &                           wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c11") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       g_momo,g_momo2,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr11(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,
     &                           g_momo,g_momo2,
     &                           wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c12") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       g_momo,g_momo2,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr12(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,
     &                           g_momo,g_momo2,
     &                           wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c13") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       g_momo,g_momo2,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr13(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,
     &                           g_momo,g_momo2,
     &                           wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1c14") then
        call wfn1_e_g_2el_mx(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                       g_momo,g_momo2,
     &                       wfn1_energy%E_2el,wfn1_energy%E_2el_hf)
        call wfn1_e_g_2el_corr14(wfn1_param,wfn1_inst,wfn1_mat,
     &                           wfn1_wave,
     &                           g_momo,g_momo2,
     &                           wfn1_energy%E_corr)
      else if (wfn1_param%corr_expr == "wfn1w") then
        call wfn1_e_g_2el_wfn(wfn1_param,wfn1_mat,wfn1_wave,
     &                        g_bfbf,g_4t1,g_4t2,
     &                        wfn1_energy%E_2el)
      else
        call errquit("invalid correlation energy option",0,UERR)
      endif
!     call ga_dgop(1400,v_dca,nmo,"+")
!     call ga_dgop(1401,v_dcb,nmo,"+")
!     call ga_dgop(1402,v_dpa,nmo,"+")
!     call ga_dgop(1403,v_dpb,nmo,"+")
c
c     The DFT correlation terms
c
      call wfn1_e_g_dft(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
!    &                  wfn1_deriv,v_dca,v_dcb,
     &                  g_bfbf,g_bfbf2,g_bfbf3,g_bfbf4,g_bfmo,g_momo,
     &                  wfn1_energy%E_dft)
cDEBUG
c     call wfn1_print_max(wfn1_deriv,"wfn1_e_g B")
cDEBUG
c
c     The 1-electron terms
c
      call wfn1_e_g_1el(wfn1_param,wfn1_mat,wfn1_wave,
!    &                  wfn1_part_drv,v_dpa,v_dpb,
     &                  g_bfmo,g_momo,wfn1_energy%E_1el)
cDEBUG
c     call wfn1_print_max(wfn1_deriv,"wfn1_e_g C")
cDEBUG
c
c     The entropy term
c
      call wfn1_e_g_ent(wfn1_param,wfn1_wave,
!    &                  v_dca,v_dcb,
     &                  wfn1_energy%E_ent)
c
      call ga_sync()
cDEBUG
c     call wfn1_print_max(wfn1_deriv,"wfn1_e_g D")
cDEBUG
c
c     Transform the correlation function Fock matrices from natural
c     orbital to correlation function basis
c
c     First add v_dca to the diagonal of the correlation function
c     Fock matrix in the natural orbital basis.
c
!     call ga_distribution(wfn1_deriv%g_dca,iproc,ilo,ihi,jlo,jhi)
!     allocate(zbuf(ilo:ihi,jlo:jhi))
!     call ga_get(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     do ii = max(ilo,jlo), min(ihi,jhi)
!       zbuf(ii,ii) = zbuf(ii,ii) + v_dca(ii)
!     enddo
!     call ga_put(wfn1_deriv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     deallocate(zbuf)
c
!     call ga_distribution(wfn1_part_drv%g_dca,iproc,ilo,ihi,jlo,jhi)
!     allocate(zbuf(ilo:ihi,jlo:jhi))
!     call ga_get(wfn1_part_drv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     do ii = max(ilo,jlo), min(ihi,jhi)
!       zbuf(ii,ii) = zbuf(ii,ii) + v_dpa(ii)
!       v_dcna(ii) = zbuf(ii,ii)
!     enddo
!     call ga_put(wfn1_part_drv%g_dca,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     deallocate(zbuf)
c
!     call ga_add(1.0d0,wfn1_deriv%g_dca,1.0d0,wfn1_part_drv%g_dca,
!    &            wfn1_deriv%g_dca)
c
c     Second transform the correlation function Fock matrix from the
c     natural orbital to the correlation function basis.
c
!     call ga_sync()
!     call ga_dgemm('T','N',nmo,nmo,nmo,1.0d0,wfn1_wave%g_ca,
!    &              wfn1_deriv%g_dca,0.0d0,g_momo)
!     call ga_dgemm('N','N',nmo,nmo,nmo,1.0d0,g_momo,wfn1_wave%g_ca,
!    &              0.0d0,wfn1_deriv%g_dca)
c
c     First add v_dcb to the diagonal of the correlation function
c     Fock matrix in the natural orbital basis.
c
!     call ga_distribution(wfn1_deriv%g_dcb,iproc,ilo,ihi,jlo,jhi)
!     allocate(zbuf(ilo:ihi,jlo:jhi))
!     call ga_get(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     do ii = max(ilo,jlo), min(ihi,jhi)
!       zbuf(ii,ii) = zbuf(ii,ii) + v_dcb(ii)
!     enddo
!     call ga_put(wfn1_deriv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     deallocate(zbuf)
c
!     call ga_distribution(wfn1_part_drv%g_dcb,iproc,ilo,ihi,jlo,jhi)
!     allocate(zbuf(ilo:ihi,jlo:jhi))
!     call ga_get(wfn1_part_drv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     do ii = max(ilo,jlo), min(ihi,jhi)
!       zbuf(ii,ii) = zbuf(ii,ii) + v_dpb(ii)
!       v_dcnb(ii) = zbuf(ii,ii)
!     enddo
!     call ga_put(wfn1_part_drv%g_dcb,ilo,ihi,jlo,jhi,zbuf,ihi-ilo+1)
!     deallocate(zbuf)
c
!     call ga_add(1.0d0,wfn1_deriv%g_dcb,1.0d0,wfn1_part_drv%g_dcb,
!    &            wfn1_deriv%g_dcb)
c
c     Second transform the correlation function Fock matrix from the
c     natural orbital to the correlation function basis.
c
!     call ga_sync()
!     call ga_dgemm('T','N',nmo,nmo,nmo,1.0d0,wfn1_wave%g_cb,
!    &              wfn1_deriv%g_dcb,0.0d0,g_momo)
!     call ga_dgemm('N','N',nmo,nmo,nmo,1.0d0,g_momo,wfn1_wave%g_cb,
!    &              0.0d0,wfn1_deriv%g_dcb)
c
c     Fetch diagonals
c
!     v_dna = 0.0d0
!     call ga_add(1.0d0,wfn1_deriv%g_dna,1.0d0,wfn1_part_drv%g_dna,
!    &            wfn1_deriv%g_dna)
!     call ga_distribution(wfn1_deriv%g_dna,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_deriv%g_dna,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dna(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
c
!     v_dnna = 0.0d0
!     call ga_distribution(wfn1_part_drv%g_dna,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_part_drv%g_dna,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dnna(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
!
!     v_dnb = 0.0d0
!     call ga_add(1.0d0,wfn1_deriv%g_dnb,1.0d0,wfn1_part_drv%g_dnb,
!    &            wfn1_deriv%g_dnb)
!     call ga_distribution(wfn1_deriv%g_dnb,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_deriv%g_dnb,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dnb(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
!
!     v_dnnb = 0.0d0
!     call ga_distribution(wfn1_part_drv%g_dnb,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_part_drv%g_dnb,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dnnb(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
!
!     v_dca = 0.0d0
!     call ga_distribution(wfn1_deriv%g_dca,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_deriv%g_dca,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dca(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
!
!     v_dcb = 0.0d0
!     call ga_distribution(wfn1_deriv%g_dcb,iproc,ilo,ihi,jlo,jhi)
!     mxlo = max(ilo,jlo)
!     mnhi = min(ihi,jhi)
!     if (mxlo.le.mnhi) then
!       allocate(buf(mxlo:mnhi,mxlo:mnhi))
!       call ga_get(wfn1_deriv%g_dcb,mxlo,mnhi,mxlo,mnhi,
!    &              buf,mnhi-mxlo+1)
!       do ii = mxlo, mnhi
!         v_dcb(ii)  = buf(ii,ii)
!       enddo
!       deallocate(buf)
!     endif
!
!     call ga_dgop(1402,v_dna, nmo,"+")
!     call ga_dgop(1403,v_dnb, nmo,"+")
!     call ga_dgop(1404,v_dca, nmo,"+")
!     call ga_dgop(1405,v_dcb, nmo,"+")
!     call ga_dgop(1406,v_dnna,nmo,"+")
!     call ga_dgop(1407,v_dnnb,nmo,"+")
!     call ga_dgop(1406,v_dcna,nmo,"+")
!     call ga_dgop(1407,v_dcnb,nmo,"+")
cDEBUG
!     call wfn1_set_entropy_eorb(v_dcna,v_dcnb,nmo)
c     call wfn1_print_orb_energies(wfn1_param,
c    &     v_dna,v_dnb,v_dca,v_dcb,
c    &     v_dnna,v_dnnb,v_dcna,v_dcnb)
cDEBUG
!     deallocate(v_dpa,v_dpb)
c
      wfn1_energy%E_tot = wfn1_energy%E_1el + wfn1_energy%E_2el
     &                  + wfn1_energy%E_nuc + wfn1_energy%E_ent
     &                  + wfn1_energy%E_dft
     &                  + wfn1_energy%E_corr
c
      call ga_sync()
c
      end subroutine wfn1_e_g
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_e_g_1el(wfn1_param,wfn1_mat,wfn1_wave,
!    &                        wfn1_deriv,v_dca,v_dcb,
     &                        g_bfmo,g_momo,E_1el)
      implicit none
C>
C> \brief Calculate the 1-electron contributions to the energy.
C>
C> The algorithm is rather simple:
C> 1. Transform the 1-electron energy matrix from the basis functions
C>    to the natural orbital basis
C> 2. Select the diagonal of the result of 1. 
C> 3. Multiply the elements from 2. with the occupation numbers and sum
C>    to get the 1-electron energy.
C> 
!include "wfn1_matcnst.fh"
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
!     type(wfn1_drv),     intent(inout) :: wfn1_deriv
c
!     double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
!     double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer, intent(in) :: g_bfmo ! nbf x nmo scratch GA
      integer, intent(in) :: g_momo ! nmo x nmo scratch GA
c
      double precision, intent(out) :: E_1el    ! the 1-electron energy
c
      double precision, allocatable :: occa(:)  ! the alpha-occupation 
                                                ! numbers
      double precision, allocatable :: occb(:)  ! the beta-occupation 
                                                ! numbers
      double precision, allocatable :: e1a(:)   ! the alpha-1-electron
                                                ! expectation energies
      double precision, allocatable :: e1b(:)   ! the beta-1-electron
                                                ! expectation energies
      double precision, allocatable :: buf(:,:) ! a local buffer
c
c     local
c
      integer iclo, ichi, jclo, jchi ! related to the correlation 
                                     ! functions
      integer imlo, imhi, jmlo, jmhi ! related to the mo x mo 
                                     ! matrix
      integer ii, jj                 ! counters
      integer nmo                    ! number of MOs
      integer nea                    ! number of alpha electrons
      integer neb                    ! number of beta  electrons
      integer iproc, nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      nmo   = wfn1_param%nmo
      nea   = wfn1_param%nea
      neb   = wfn1_param%neb
c
      allocate(occa(1:nmo),occb(1:nmo))
      allocate(e1a(1:nmo),e1b(1:nmo))
c
      occa = 0.0d0
      occb = 0.0d0
      e1a  = 0.0d0
      e1b  = 0.0d0
c
c     Do the ALPHA spin channel
c     -------------------------
c
c     Calculate the alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
c
c     Do step 1.
c
      call ga_dgemm('n','n',wfn1_param%nbf,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_mat%g_tv,wfn1_wave%g_na,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',wfn1_param%nmo,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_wave%g_na,g_bfmo,
     &              0.0d0,g_momo)
c
c     Do step 2.
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1a(ii),1)
      enddo
      call ga_dgop(1406,e1a,wfn1_param%nmo,"+")
c
c     Do step 3.
c
      do ii = 1, wfn1_param%nmo
        E_1el = E_1el + occa(ii)*e1a(ii)
      enddo
c
c     Do step 4.
c
!     allocate(buf(imlo:imhi,jmlo:jmhi))
!     call ga_get(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
!     do jj = jmlo, jmhi
!       do ii = imlo, imhi
!         buf(ii,jj) = buf(ii,jj)*occa(ii)
!       enddo
!     enddo
!     call ga_put(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
!     deallocate(buf)
!     call ga_add(1.0d0,wfn1_deriv%g_dna,
!    &            1.0d0,g_momo,wfn1_deriv%g_dna)
c
c     Do step 5.
c
!     v_dca = v_dca + e1a
c
c     Do the BETA spin channel
c     ------------------------
c
c     Calculate the beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
c
c     Do step 1.
c
      call ga_dgemm('n','n',wfn1_param%nbf,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_mat%g_tv,wfn1_wave%g_nb,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',wfn1_param%nmo,wfn1_param%nmo,
     &              wfn1_param%nbf,1.0d0,wfn1_wave%g_nb,g_bfmo,
     &              0.0d0,g_momo)
c
c     Do step 2.
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1b(ii),1)
      enddo
      call ga_dgop(1407,e1b,wfn1_param%nmo,"+")
c
c     Do step 3.
c
      do ii = 1, wfn1_param%nmo
        E_1el = E_1el + occb(ii)*e1b(ii)
      enddo
c
c     Do step 4.
c
!     allocate(buf(imlo:imhi,jmlo:jmhi))
!     call ga_get(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
!     do jj = jmlo, jmhi
!       do ii = imlo, imhi
!         buf(ii,jj) = buf(ii,jj)*occb(ii)
!       enddo
!     enddo
!     call ga_put(g_momo,imlo,imhi,jmlo,jmhi,buf,imhi-imlo+1)
!     deallocate(buf)
!     call ga_add(1.0d0,wfn1_deriv%g_dnb,
!    &            1.0d0,g_momo,wfn1_deriv%g_dnb)
c
c     Do step 5.
c
!     v_dcb = v_dcb + e1b
c
c     Done the work, now tidy up
c
      deallocate(e1a,e1b)
      deallocate(occa,occb)
      call ga_sync()
      end
c
c-----------------------------------------------------------------------
c
!>
!> \brief Calculate the entropy derived energy contribution
!>
!> The entropy derived contribution is given by
!> \f{eqnarray}{
!>   E_{\mathrm{ent}} &=& -T_{\mathrm{bath}}S \\\\
!>   S &=& -\sum_i \rho_i\ln(\rho_i) + (1-\rho_i)ln(1-\rho_i)
!> \f}
!> where \f$T_{\mathrm{bath}}\f$ is the bath temperature, \f$S\f$ is
!> the entropy, and \f$\rho_i\f$ is the occupation number of the
!> \f$i\f$-th natural orbital.
!>
!> The derivative of the entropy with respect to the occupation numbers
!> is given by
!> \f{eqnarray}{
!>   \frac{\partial S}{\partial\rho_i}
!>   &=& -\ln(\rho_i)+ln(1-\rho_i)
!> \f}
!> Obviously if \f$\rho_i = 1\f$ or \f$\rho_i = 0\f$ both the entropy
!> and its derivative are undefined. Hence we need to introduce a 
!> cut-off \f$\rho_c\f$. 
!>
!> If \f$\rho_i \lt \rho_c\f$ then 
!> \f{eqnarray}{
!>   S(\rho_i) &=& 0 \\\\
!>   \frac{\partial S}{\partial\rho_i} 
!>   &=& -\ln(\rho_c)
!> \f}
!> If \f$\rho_i \gt 1-\rho_c\f$ then 
!> \f{eqnarray}{
!>   S(\rho_i) &=& 0 \\\\
!>   \frac{\partial S}{\partial\rho_i} 
!>   &=& \ln(1-\rho_c)
!> \f}
!>
      subroutine wfn1_e_g_ent(wfn1_param,wfn1_wave,
!    &                        v_dca,v_dcb,
     &                        E_ent)
      implicit none
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
#include "global.fh"
#include "mafdecls.fh"
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
c
!     double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
!     double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
      double precision, intent(inout) :: E_ent
c
      double precision, allocatable :: occa(:)  ! the alpha-occupation
                                                ! numbers
      double precision, allocatable :: occb(:)  ! the beta-occupation
                                                ! numbers
      double precision, allocatable :: e1a(:)   ! the alpha-1-electron
                                                ! energies
      double precision, allocatable :: e1b(:)   ! the beta-1-electron
                                                ! energies
      double precision, allocatable :: buf(:,:) ! a local buffer
c
c     local
c
      double precision :: dr, d1r    ! entropy derivatives
      double precision :: rho_c      ! cut-off for occupation numbers
c
      integer iclo, ichi, jclo, jchi ! related to the correlation
                                     ! functions
      integer ii, jj                 ! counters
      integer iproc                  ! the current processor rank
      integer nproc                  ! the total number of processor
                                     ! ranks
      integer :: ncorea ! the number of alpha core states
      integer :: ncoreb ! the number of beta  core states
      integer :: nvirta ! the number of alpha virtual states
      integer :: nvirtb ! the number of beta  virtual states
      integer :: nea    ! the number of alpha electrons
      integer :: neb    ! the number of beta  electrons
      integer :: nmo    ! the number of molecular orbitals
c
      double precision :: s  ! The entropy of a natural orbital
      double precision :: ds ! The derivative
c
      if (wfn1_param%Tbath.le.0.0d0) return
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      allocate(occa(1:wfn1_param%nmo),occb(1:wfn1_param%nmo))
!     allocate(e1a(1:wfn1_param%nmo),e1b(1:wfn1_param%nmo))
c
      ncorea = wfn1_param%ncorea
      ncoreb = wfn1_param%ncoreb
      nvirta = wfn1_param%nvirta
      nvirtb = wfn1_param%nvirtb
      nea    = wfn1_param%nea
      neb    = wfn1_param%neb
      nmo    = wfn1_param%nmo
c
      rho_c = 1.0d-10
      occa  = 0.0d0
      occb  = 0.0d0
!     e1a   = 0.0d0
!     e1b   = 0.0d0
c
c     Do the ALPHA spin channel
c     -------------------------
c
c     Calculate the alpha-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
c
c     Calculate the entropy
c
      do ii = ncorea+1, nmo-nvirta
        call wfn1_entropy1(occa(ii),s,ds)
        E_ent   = E_ent   + wfn1_param%Tbath*s
!       e1a(ii) = e1a(ii) + wfn1_param%Tbath*ds
      enddo
c
c     Do the BETA spin channel
c     ------------------------
c
c     Calculate the beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
c
c     Calculate the entropy
c
      do ii = ncoreb+1, nmo-nvirtb
        call wfn1_entropy1(occb(ii),s,ds)
        E_ent   = E_ent   + wfn1_param%Tbath*s
!       e1b(ii) = e1b(ii) + wfn1_param%Tbath*ds
      enddo
c
c     Store the results
c     -----------------
c
!     v_dca = v_dca + e1a
!     v_dcb = v_dcb + e1b
c
c     Clean up memory
c
      deallocate(occa,occb)
!     deallocate(e1a,e1b)
c
      end
c
c-----------------------------------------------------------------------
c
  
