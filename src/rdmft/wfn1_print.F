c-----------------------------------------------------------------------
c
      subroutine wfn1_print_wave(wfn1_wave)
      implicit none
c
!include "wfn1_wfn.fh"
#include "global.fh"
#include "stdio.fh"
c
      type(wfn1_wfn), intent(in) :: wfn1_wave ! The wavefunction
c
      integer iproc ! the process rank
      integer nproc ! the number of process ranks
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
        write(LuOut,'(" The alpha natural orbitals")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_na)
      if (iproc.eq.0) then
        call util_flush(LuOut)
        write(LuOut,*)
        write(LuOut,'(" The beta natural orbitals")')
        write(LuOut,'(" -------------------------")')
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_nb)
      if (iproc.eq.0) then
        call util_flush(LuOut)
        write(LuOut,*)
        write(LuOut,'(" The alpha correlation functions")')
        write(LuOut,'(" -------------------------------")')
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_ca)
      if (iproc.eq.0) then
        call util_flush(LuOut)
        write(LuOut,*)
        write(LuOut,'(" The beta correlation functions")')
        write(LuOut,'(" ------------------------------")')
        write(LuOut,*)
        call util_flush(LuOut)
      endif
      call ga_print(wfn1_wave%g_cb)
      if (iproc.eq.0) then
        call util_flush(LuOut)
      endif
c
      end subroutine wfn1_print_wave
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_deriv(wfn1_deriv)
      implicit none
c
!include "wfn1_drv.fh"
#include "global.fh"
#include "stdio.fh"
c
      type(wfn1_drv), intent(in) :: wfn1_deriv ! The derivative
c
      integer iproc ! the process rank
      integer nproc ! the number of process ranks
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      if (iproc.eq.0) then
        write(LuOut,'(" The alpha natural orbital derivatives")')
        write(LuOut,'(" -------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dna)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta natural orbitals derivatives")')
        write(LuOut,'(" -------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dnb)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The alpha correlation function derivatives")')
        write(LuOut,'(" ------------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dca)
      if (iproc.eq.0) then
        write(LuOut,*)
        write(LuOut,'(" The beta correlation function derivatives")')
        write(LuOut,'(" -----------------------------------------")')
        write(LuOut,*)
      endif
      call ga_print(wfn1_deriv%g_dcb)
c
      end subroutine wfn1_print_deriv
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_energy(wfn1_energy)
      implicit none
!include "wfn1_en.fh"
#include "global.fh"
#include "stdio.fh"
      type(wfn1_en), intent(in) :: wfn1_energy
      integer iproc
      integer nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      if (iproc.eq.0) then
        write(LuOut,'(" WFN1                     :")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,'(" Total WFN1 energy  = ",f22.10)')
     &              wfn1_energy%E_tot
        write(LuOut,'(" 1-Electron energy  = ",f22.10)')
     &              wfn1_energy%E_1el
        write(LuOut,'(" 2-Electron energy  = ",f22.10)')
     &              wfn1_energy%E_2el
        write(LuOut,'(" Correlation energy = ",f22.10)')
     &              wfn1_energy%E_corr
        write(LuOut,'(" DFT XC energy      = ",f22.10)')
     &              wfn1_energy%E_dft
        write(LuOut,'(" El. redistribution = ",f22.10)')
     &              wfn1_energy%E_redist
        write(LuOut,'(" Entropy energy     = ",f22.10)')
     &              wfn1_energy%E_ent
        write(LuOut,'(" Nuclear energy     = ",f22.10)')
     &              wfn1_energy%E_nuc
        write(LuOut,*)
        write(LuOut,'(" conventional Hartree-Fock:")')
        write(LuOut,'(" --------------------------")')
        write(LuOut,'(" Total SCF energy  = ",f22.10)')
     &        wfn1_energy%E_tot-wfn1_energy%E_2el+wfn1_energy%E_2el_hf
        write(LuOut,'(" 2-Electron energy = ",f22.10)')
     &        wfn1_energy%E_2el_hf
      endif
c
      end subroutine wfn1_print_energy
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_error(wfn1_err)
      implicit none
!include "wfn1_error.fh"
#include "global.fh"
#include "stdio.fh"
      type(wfn1_error), intent(in) :: wfn1_err
      integer :: iproc
      integer :: nproc
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      if (iproc.eq.0) then
        write(LuOut,'(" Tot:",f16.6," Na:",e16.6," Nb:",e16.6,
     &                              " Ca:",e16.6," Cb:",e16.6)')
     &    wfn1_err%total,wfn1_err%errna,wfn1_err%errnb,wfn1_err%errca,
     &    wfn1_err%errcb
      endif
c
      end subroutine wfn1_print_error
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_mat(dmat,nr,nc)
      implicit none
c
c     A simple serial matrix print routine
c
#include "global.fh"
#include "stdio.fh"
c
      integer, intent(in) :: nr ! the number of rows
      integer, intent(in) :: nc ! the number of columns
c
      double precision, intent(in) :: dmat(nr,nc)
c
      integer :: ir, ic, kk ! counters
      integer :: iproc, nproc
c
      integer, parameter :: mxc = 8 ! the maximum number of columns
                                    ! printed at a time
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      if (iproc.eq.0) then
c
        do ic = 1, nc, mxc
          write(LuOut,*)
          write(LuOut,'(5x,$)')
          do kk = ic, min(ic+mxc-1,nc)
            write(LuOut,'(i12,$)')kk
          enddo
          write(LuOut,*)
          write(LuOut,*)
          do ir = 1, nr
            write(LuOut,'(i4,$)')ir
            do kk = ic, min(ic+mxc-1,nc)
              write(LuOut,'(f12.6,$)')dmat(ir,kk)
            enddo
            write(LuOut,*)
          enddo
        enddo
c
      endif
c
      end subroutine wfn1_print_mat
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_mat2(dmat,nbf)
      implicit none
c
c     A simple serial 2RDM print routine
c
#include "global.fh"
#include "stdio.fh"
c
      integer, intent(in) :: nbf ! the number of basis functions
c
      double precision, intent(in) :: dmat(nbf,nbf,nbf,nbf)
c
      integer :: ia, ib, ic, id, kk, ll, k2 ! counters
      integer :: iproc, nproc
c
      integer, parameter :: mxc = 8 ! the maximum number of columns
                                    ! printed at a time
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      if (iproc.eq.0) then
c
        do kk = 1, nbf*nbf, mxc
          write(LuOut,*)
          write(LuOut,'("              ",$)')
          do k2 = kk, min(kk+mxc-1,nbf*nbf)
            id = mod(k2-1,nbf)+1
            ic = (k2+nbf-1)/nbf
            write(LuOut,'(" (",i3,",",i3,")    ",$)')ic,id
          enddo
          write(LuOut,10)
10        format(" (r'1,r'2)")
          write(LuOut,*)
          do ll = 1, nbf*nbf
            ib = mod(ll-1,nbf)+1
            ia = (ll+nbf-1)/nbf
            write(LuOut,'(" (",i3,",",i3,")   ",$)')ia,ib
            do k2 = kk, min(kk+mxc-1,nbf*nbf)
              id = mod(k2-1,nbf)+1
              ic = (k2+nbf-1)/nbf
              write(LuOut,'(e14.6,$)')dmat(ia,ic,ib,id)
            enddo
            write(LuOut,*)
          enddo
          write(LuOut,'(" ( r1, r2)")')
        enddo
c
      endif
c
      end subroutine wfn1_print_mat2
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_print_orb_energies(wfn1_param,e_na,e_nb,e_ca,e_cb,
     &           e_nna,e_nnb,e_cna,e_cnb)
      implicit none
#include "stdio.fh"
!include "wfn1_param.fh"
!>    Calculation parameters 
      type(wfn1_prm), intent(in)   :: wfn1_param
!>    \f$\alpha\f$ natural orbital Fock matrix eigenvalues
      double precision, intent(in) :: e_na(1:wfn1_param%nmo)
!>    \f$\beta\f$ natural orbital Fock matrix eigenvalues
      double precision, intent(in) :: e_nb(1:wfn1_param%nmo)
!>    \f$\alpha\f$ natural orbital expectation values of the 
!>    correlation function Fock matrix
      double precision, intent(in) :: e_ca(1:wfn1_param%nmo)
!>    \f$\alpha\f$ natural orbital expectation values of the 
!>    correlation function Fock matrix
      double precision, intent(in) :: e_cb(1:wfn1_param%nmo)
!>    \f$\alpha\f$ natural orbital Fock matrix eigenvalues
      double precision, intent(in) :: e_nna(1:wfn1_param%nmo)
!>    \f$\beta\f$ natural orbital Fock matrix eigenvalues
      double precision, intent(in) :: e_nnb(1:wfn1_param%nmo)
!>    \f$\alpha\f$ natural orbital expectation values of the 
!>    correlation function Fock matrix
      double precision, intent(in) :: e_cna(1:wfn1_param%nmo)
!>    \f$\alpha\f$ natural orbital expectation values of the 
!>    correlation function Fock matrix
      double precision, intent(in) :: e_cnb(1:wfn1_param%nmo)
c
      double precision :: e1_na(1:wfn1_param%nmo)
      double precision :: e1_nb(1:wfn1_param%nmo)
      double precision :: e1_ca(1:wfn1_param%nmo)
      double precision :: e1_cb(1:wfn1_param%nmo)
      double precision :: e1_nna(1:wfn1_param%nmo)
      double precision :: e1_nnb(1:wfn1_param%nmo)
      double precision :: e1_cna(1:wfn1_param%nmo)
      double precision :: e1_cnb(1:wfn1_param%nmo)
      integer :: ii
c
      e1_na=e_na
      e1_nb=e_nb
      e1_ca=e_ca
      e1_cb=e_cb
      e1_nna=e_nna
      e1_nnb=e_nnb
      e1_cna=e_cna
      e1_cnb=e_cnb
      call wfn1_sort(e1_na,wfn1_param%nmo)
      call wfn1_sort(e1_nb,wfn1_param%nmo)
      call wfn1_sort(e1_ca,wfn1_param%nmo)
      call wfn1_sort(e1_cb,wfn1_param%nmo)
c     call wfn1_sort(e1_nna,wfn1_param%nmo)
c     call wfn1_sort(e1_nnb,wfn1_param%nmo)
c     call wfn1_sort(e1_cna,wfn1_param%nmo)
c     call wfn1_sort(e1_cnb,wfn1_param%nmo)
      write(LuOut,'(12x,"  Alpha ",4x,"  Beta  ",8x,
     &                  "  Alpha ",4x,"  Beta  ")')
      write(LuOut,'(12x,"--------",4x,"--------",8x,
     &                  "--------",4x,"--------")')
      do ii = 1, wfn1_param%nmo
        write(LuOut,'(i8,2f12.4,4x,2f12.4)')ii,e1_na(ii),e1_nb(ii),
     &                                      e1_ca(ii),e1_cb(ii)
      enddo
cDEBUG
      do ii = 1, wfn1_param%nmo
        write(75,'(2x,f24.16,$)')e1_na(ii)
        write(76,'(2x,f24.16,$)')e1_ca(ii)
        write(77,'(2x,f24.16,$)')e1_nna(ii)
        write(78,'(2x,f24.16,$)')e1_cna(ii)
        write(89,'(2x,f24.16,$)')e1_na(ii)-e1_nna(ii)
        write(90,'(2x,f24.16,$)')e1_ca(ii)-e1_cna(ii)
      enddo
      write(75,*)
      write(76,*)
      write(77,*)
      write(78,*)
      write(89,*)
      write(90,*)
      close(75)
      close(76)
      close(77)
      close(78)
      close(89)
      close(90)
cDEBUG
      end subroutine wfn1_print_orb_energies
c
c-----------------------------------------------------------------------
c
!> \brief Sort the values in A in increasing order
!>
      subroutine wfn1_sort(a,n)
      implicit none
      integer, intent(in) :: n
      double precision, intent(inout) :: a(n)
      integer :: ii, jj, itmp
      double precision tmp
      do ii = 1, n
        tmp = a(ii)
        itmp = ii
        do jj = ii+1,n
          if (a(jj).lt.tmp) then
            tmp = a(jj)
            itmp = jj
          endif
        enddo
        if (itmp.ne.ii) then
          a(itmp) = a(ii)
          a(ii) = tmp
        endif
      enddo
      end subroutine wfn1_sort
c
c-----------------------------------------------------------------------
