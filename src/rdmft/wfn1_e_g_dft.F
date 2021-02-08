C> \brief Calculate the DFT contributions to the energy,
C> the natural orbital and the correlation function Fock matrices
C>
C> The 1-electron term to the Fock matrices will be calculated as:
C> - The natural orbital term will calculated in the natural orbital
C>   basis
C> - The correlation function Fock matrices will be calculated in the
C>   natural orbital basis because the matrix is diagonal in that 
C>   basis
C> Regarding the correlation function Fock matrices these will be
C> calculated in the natural orbital basis in this routine. Once the
C> full diagonal elements have calculated (including the 2-electron
C> terms) then the matrix will be transformed to the correlation
C> function basis. Doing this transformation prematurely just
C> generates unnecesary work.
C> The algorithm is rather simple:
C> 1. Transform the 1-electron energy matrix from the basis functions
C>    to the natural orbital basis
C> 2. Select the diagonal of the result of 1. 
C> 3. Multiply the elements from 2. with the occupation numbers and sum
C>    to get the 1-electron energy.
C> 4. Scale the result from 1. with symmetrized occupation number 
C>    factors to get the natural orbital Fock matrix 1-electron 
C>    contribution term.
C> 5. Add the data from 2. to the diagonal of the correlation function
C>    Fock matrix in the natural orbital basis.
C> 
      subroutine wfn1_e_g_dft(wfn1_param,wfn1_inst,wfn1_mat,wfn1_wave,
     &                        wfn1_deriv,
     &                        v_dca,v_dcb,
     &                        g_bfbf_a,g_bfbf_b,g_bfbf_xa,g_bfbf_xb,
     &                        g_bfmo,g_momo,E_dft)
      implicit none
!include "wfn1_matcnst.fh"
!include "wfn1_param.fh"
!include "wfn1_wfn.fh"
!include "wfn1_drv.fh"
#include "global.fh"
#include "mafdecls.fh"
      target :: dbl_mb
c
      type(wfn1_prm),     intent(in)    :: wfn1_param
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
      type(wfn1_drv),     intent(inout) :: wfn1_deriv
c
      double precision, intent(inout) :: v_dca(1:wfn1_param%nmo)
      double precision, intent(inout) :: v_dcb(1:wfn1_param%nmo)
c
      integer, intent(in) :: g_bfbf_a  ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf_b  ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf_xa ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfbf_xb ! nbf x nbf scratch GA
      integer, intent(in) :: g_bfmo    ! nbf x nmo scratch GA
      integer, intent(in) :: g_momo    ! nmo x nmo scratch GA
c
      double precision, intent(out) :: E_dft    ! the 1-electron energy
c
      double precision, allocatable :: occa(:)  ! the alpha-occupation 
                                                ! numbers
      double precision, allocatable :: occb(:)  ! the beta-occupation 
                                                ! numbers
      double precision, allocatable :: e1a(:)   ! the alpha-1-electron
                                                ! energies
      double precision, allocatable :: e1b(:)   ! the beta-1-electron
                                                ! energies
c
      double precision, allocatable :: buf (:,:) ! a local buffer
      double precision, allocatable :: bufi(:,:) ! a local buffer
      double precision, allocatable :: bufj(:,:) ! a local buffer
c
      double precision, pointer :: dd(:,:) ! point for the density
                                           ! matrix
c
c     functions
c
      logical, external :: xc_gotxc ! has a density functional been 
                                    ! specified?
c
c     local
c
      double precision :: rho_n      ! total number of electrons
      double precision :: exc, dume, dum1, dum2 ! junk
      integer iclo, ichi, jclo, jchi ! related to the correlation 
                                     ! functions
      integer imlo, imhi, jmlo, jmhi ! related to the mo x mo 
                                     ! matrix
      integer iblo, ibhi, jblo, jbhi ! related to the bf x bf 
                                     ! matrix
      integer idx, ldx
      integer na, nb, nn             ! dimensions
      integer ii, jj                 ! counters
      integer nbf                    ! number of basis functions
      integer nmo                    ! number of MOs
      integer nea                    ! number of alpha electrons
      integer neb                    ! number of beta  electrons
      integer iproc, nproc
c
      integer :: g_dens(2)
      integer :: g_vxc(2)
c
      E_dft = 0.0d0
cDEBUG
c     write(*,*)'HvD: xc_gotxc = ',xc_gotxc()
cDEBUG
      if (.not.xc_gotxc()) return
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      nbf   = wfn1_param%nbf
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
      call ga_zero(g_bfbf_a)
      call ga_zero(g_bfbf_b)
c
c     Calculate the alpha-occupation numbers
c
cDEBUG
c     write(*,*)'HVD: nea   = ',nea
c     write(*,*)'HVD: neb   = ',neb
c     write(*,*)'HVD: icnfa = ',wfn1_wave%icnfa
c     write(*,*)'HVD: icnfb = ',wfn1_wave%icnfb
c     write(*,*)'HVD: dcnfa = ',wfn1_wave%dcnta
c     write(*,*)'HVD: dcnfb = ',wfn1_wave%dcntb
c     call util_flush(6)
c     call ga_print(wfn1_wave%g_ca)
c     call ga_print(wfn1_wave%g_cb)
c     call util_flush(6)
cDEBUG
      call wfn1_calc_occ(wfn1_wave%g_ca,nmo,nea,wfn1_wave%icnfa,
     &                   wfn1_wave%dcnta,occa)
c
c     Calculate the beta-occupation numbers
c
      call wfn1_calc_occ(wfn1_wave%g_cb,nmo,neb,wfn1_wave%icnfb,
     &                   wfn1_wave%dcntb,occb)
c
c     Calculate the alpha-density matrix
c
      call ga_distribution(g_bfbf_a,iproc,iblo,ibhi,jblo,jbhi)
      call ga_access(g_bfbf_a,iblo,ibhi,jblo,jbhi,idx,ldx)
      na = (ibhi-iblo+1)
      nb = (jbhi-jblo+1)
      nn = na*nb-1
      dd(iblo:ibhi,jblo:jbhi) => dbl_mb(idx:idx+nn)
      allocate(bufi(iblo:ibhi,1:nmo))
      allocate(bufj(jblo:jbhi,1:nmo))
      call ga_get(wfn1_wave%g_na,iblo,ibhi,1,nmo,bufi,na)
      call ga_get(wfn1_wave%g_na,jblo,jbhi,1,nmo,bufj,nb)
      do jj = jblo, jbhi
        do ii = iblo, ibhi
          do nn = 1, nmo
            dd(ii,jj) = dd(ii,jj) + bufi(ii,nn)*bufj(jj,nn)*occa(nn)
          enddo
        enddo
      enddo
      deallocate(bufi,bufj)
c
c     Calculate the beta-density matrix
c
      call ga_distribution(g_bfbf_b,iproc,iblo,ibhi,jblo,jbhi)
      call ga_access(g_bfbf_b,iblo,ibhi,jblo,jbhi,idx,ldx)
      na = (ibhi-iblo+1)
      nb = (jbhi-jblo+1)
      nn = na*nb-1
      dd(iblo:ibhi,jblo:jbhi) => dbl_mb(idx:idx+nn)
      allocate(bufi(iblo:ibhi,1:nmo))
      allocate(bufj(jblo:jbhi,1:nmo))
      call ga_get(wfn1_wave%g_nb,iblo,ibhi,1,nmo,bufi,na)
      call ga_get(wfn1_wave%g_nb,jblo,jbhi,1,nmo,bufj,nb)
      do jj = jblo, jbhi
        do ii = iblo, ibhi
          do nn = 1, nmo
            dd(ii,jj) = dd(ii,jj) + bufi(ii,nn)*bufj(jj,nn)*occb(nn)
          enddo
        enddo
      enddo
      deallocate(bufi,bufj)
cDEBUG
c     write(*,*)'HVD: occa = ',occa
c     write(*,*)'HVD: occb = ',occb
c     call util_flush(6)
c     call ga_print(g_bfbf_a)
c     call ga_print(g_bfbf_b)
c     call util_flush(6)
cDEBUG
c
c     Call the DFT module
c
      g_dens(1) = g_bfbf_a
      g_dens(2) = g_bfbf_b
      g_vxc(1)  = g_bfbf_xa
      g_vxc(2)  = g_bfbf_xb
      call xc_getv(wfn1_inst%rtdb,exc,dume,1,0,0,g_dens,g_vxc,
     &             .false.,0,0,0,.false.,rho_n,dum1,dum2,0)
c     if (abs(rho_n-dble(nea+neb)).gt.1.0d-4) then
c       write(*,*)'No electrons problem: ',nea+neb,rho_n
c     endif
      E_dft = ga_ddot(g_dens(1),g_vxc(1))
     &      + ga_ddot(g_dens(2),g_vxc(2))
c
c     Generate proper alpha-Fock matrix contributions
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,g_bfbf_xa,wfn1_wave%g_na,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_na,g_bfmo,
     &              0.0d0,g_momo)
c
c     Extract the diagonal elements for the alpha correlation function Fock
c     matrix
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1a(ii),1)
      enddo
      call ga_dgop(1406,e1a,nmo,"+")
c
c     Scale Fock matrix elements with occupation numbers for the
c     alpha natural orbital Fock matrix
c
      call ga_access(g_momo,imlo,imhi,jmlo,jmhi,idx,ldx)
      nn = (imhi-imlo-1)*(jmhi-jmlo-1)-1
      dd(imlo:imhi,jmlo:jmhi) => dbl_mb(idx:idx+nn)
      do jj = jmlo, jmhi
        do ii = imlo, imhi
          dd(ii,jj) = dd(ii,jj)*occa(ii)
        enddo
      enddo
      call ga_sync()
      call ga_add(1.0d0,wfn1_deriv%g_dna,
     &            1.0d0,g_momo,wfn1_deriv%g_dna)
c
c     Add the alpha correlation function Fock matrix contributions
c
      v_dca = v_dca + e1a
c
c     Generate proper beta-Fock matrix contributions
c
      call ga_dgemm('n','n',nbf,nmo,nbf,1.0d0,g_bfbf_xb,wfn1_wave%g_nb,
     &              0.0d0,g_bfmo)
      call ga_dgemm('t','n',nmo,nmo,nbf,1.0d0,wfn1_wave%g_nb,g_bfmo,
     &              0.0d0,g_momo)
c
c     Extract the diagonal elements for the beta correlation function Fock
c     matrix
c
      call ga_distribution(g_momo,iproc,imlo,imhi,jmlo,jmhi)
      do ii = max(imlo,jmlo), min(imhi,jmhi)
        call ga_get(g_momo,ii,ii,ii,ii,e1b(ii),1)
      enddo
      call ga_dgop(1407,e1b,nmo,"+")
c
c     Scale Fock matrix elements with occupation numbers for the
c     beta natural orbital Fock matrix
c
      call ga_access(g_momo,imlo,imhi,jmlo,jmhi,idx,ldx)
      nn = (imhi-imlo-1)*(jmhi-jmlo-1)-1
      dd(imlo:imhi,jmlo:jmhi) => dbl_mb(idx:idx+nn)
      do jj = jmlo, jmhi
        do ii = imlo, imhi
          dd(ii,jj) = dd(ii,jj)*occb(ii)
        enddo
      enddo
      call ga_sync()
      call ga_add(1.0d0,wfn1_deriv%g_dnb,
     &            1.0d0,g_momo,wfn1_deriv%g_dnb)
c
c     Add the beta correlation function Fock matrix contributions
c
      v_dcb = v_dcb + e1b
c
c     Clean up
c
      deallocate(e1a,e1b)
      deallocate(occa,occb)
      call ga_sync()
      end
