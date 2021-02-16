c-----------------------------------------------------------------------
C>
C> \brief Calculate occupation numbers from a set of correlation
C> functions
C>
C> Note that now the correlation functions are complex valued
C> but the occupation numbers are still real because we are
C> taking the norm.
C>
      subroutine wfn1_calc_occ(g_c,nmo,ne,icnf,dcnt,occ)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer, intent(in) :: g_c      ! the correlation functions
      integer, intent(in) :: nmo      ! the number of molecular orbitals
      integer, intent(in) :: ne       ! the number of electrons
      integer, intent(in) :: icnf(nmo)! the electron configuration
c
      double precision, intent(in)  :: dcnt(nmo)! the electron fraction
c
      double precision, intent(out) :: occ(nmo) ! the occupation numbers
c
      integer iproc ! the processor rank
      integer nproc ! the number of processor ranks
c
      integer ilo, ihi, jlo, jhi ! patch limits
      integer ii, ir             ! counters
      integer ityp               ! type of g_c
c
      complex(kind=8), allocatable :: buf(:,:)
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      call ga_inquire(g_c,ityp,ii,ir)
      if (ityp.ne.MT_DCPL) then
        call errquit("wfn1_calc_occ: g_c not complex",0,UERR)
      endif
c
      occ = 0.0d0
      call ga_distribution(g_c,iproc,ilo,ihi,jlo,jhi)
      allocate(buf(ilo:ihi,jlo:jhi))
      if (ihi.ge.ilo.and.jhi.ge.jlo) then
        call ga_get(g_c,ilo,ihi,jlo,jhi,buf,ihi-ilo+1)
      else
c       This process does not hold any part of g_c
      endif
      do ir = 1, ne
        if(icnf(ir).ge.jlo.and.icnf(ir).le.jhi)
     &  then
          do ii = ilo, ihi
            occ(ii) = occ(ii)
     &              + buf(ii,icnf(ir)) * conjg(buf(ii,icnf(ir))) *
     &                dcnt(ir)
          enddo
        endif
      enddo
      deallocate(buf)
      call ga_dgop(1200,occ,nmo,"+")
c
      end subroutine wfn1_calc_occ
c
c-----------------------------------------------------------------------
C>
C> \brief Calculate the 2-electron occupation numbers from a set of
C> correlation functions
C>
C> Note the correlation functions are complex valued now. The occupation
C> numbers may remain real because we are taking norms.
C>
!     subroutine wfn1_calc_occ2(g_c,nmo,ne,icnf,dcnt,g_occ)
!     implicit none
!include "global.fh"
!include "mafdecls.fh"
!include "errquit.fh"
c
!     integer, intent(in) :: g_c       ! the correlation functions
!     integer, intent(in) :: nmo       ! the number of molecular orbitals
!     integer, intent(in) :: ne        ! the number of electrons
!     integer, intent(in) :: icnf(nmo) ! the electron configuration
!     integer, intent(in) :: g_occ     ! the 2-electron occupation 
!                                      ! numbers
c
!     double precision, intent(in) :: dcnt(nmo)! the electron fraction
c
!     integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
!     complex(kind=8), allocatable :: bocc(:,:)  ! buffer for 
!                                                ! occupation numbers
!     complex(kind=8), allocatable :: bcori(:,:) ! buffer for 
!                                                ! correlation function
!     complex(kind=8), allocatable :: bcorj(:,:) ! buffer for 
!                                                ! correlation function
c
!     integer iproc ! the current processor rank
!     integer nproc ! the number of processor ranks
c
!     integer is, it ! counters for the correlation functions
!     integer ii, ij ! counters for the natural orbitals
c
!     nproc = ga_nnodes()
!     iproc = ga_nodeid()
c
!     call ga_distribution(g_occ,iproc,iolo,iohi,jolo,johi)
!     allocate(bocc(iolo:iohi,jolo:johi))
!     allocate(bcori(iolo:iohi,1:nmo))
!     allocate(bcorj(jolo:johi,1:nmo))
!     bocc = 0.0d0
!     call ga_get(g_c,iolo,iohi,1,nmo,bcori,iohi-iolo+1)
!     call ga_get(g_c,jolo,johi,1,nmo,bcorj,johi-jolo+1)
!     do ij = jolo, johi
!       do ii = iolo, iohi
!         do it = 1, ne
!           do is = 1, ne
!             bocc(ii,ij) = bocc(ii,ij) 
!    +                    + bcori(ii,icnf(is))*conjg(bcori(ii,icnf(is)))
!    +                    * bcorj(ij,icnf(it))*conjg(bcorj(ij,icnf(it)))
!    +                    * dcnt(is)*dcnt(it)
!    +                    - bcori(ii,icnf(is))*conjg(bcori(ii,icnf(it)))
!    +                    * bcorj(ij,icnf(it))*conjg(bcorj(ij,icnf(is)))
!    +                    * dcnt(is)*dcnt(it)
!           enddo
!         enddo
!       enddo
!     enddo
!     call ga_put(g_occ,iolo,iohi,jolo,johi,bocc,iohi-iolo+1)
c
!     deallocate(bcorj)
!     deallocate(bcori)
!     deallocate(bocc)
c
!     end subroutine wfn1_calc_occ2
c
c-----------------------------------------------------------------------
c
C>
C> \brief Calculate the 2-electron exchange occupation numbers
C>
C> These quantities are the
C> \f$-\sum_{s,t}C_{is}C_{it}^*C_{jt}C_{js}^*\f$ terms in the 2-electron
C> occupation numbers.
C>
C> Note that for the 2-electron density matrix the Coulomb and
C> exchange terms are separated. The reason is that the correlation
C> function Fock matrices for the Coulomb and exchange terms have to be
C> handled differently because of their nature.
C>
C> Note that although the coefficients \f$C_{is}\f$ are complex valued
C> one should expect the occupation numbers to remain real as for every
C> term its complex conjugate will also be added. Hence the imaginary
C> parts sum up to zero.
C>
      subroutine wfn1_calc_ocx(g_c,nmo,ne,icnf,dcnt,g_ocx)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "errquit.fh"
c
      integer, intent(in) :: g_c       ! the correlation functions
      integer, intent(in) :: nmo       ! the number of molecular orbitals
      integer, intent(in) :: ne        ! the number of electrons
      integer, intent(in) :: icnf(nmo) ! the electron configuration
c
      double precision, intent(in) :: dcnt(nmo) ! the electron fraction
c
      integer, intent(in) :: g_ocx  ! the 2-electron exchange 
                                    ! occupation numbers
c
      integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
      real(kind=8),    allocatable :: bocc(:,:)  ! buffer for 
                                                 ! occupation numbers
      complex(kind=8), allocatable :: boccz(:,:) ! buffer for 
                                                 ! occupation numbers
      complex(kind=8), allocatable :: bcori(:,:) ! buffer for 
                                                 ! correlation function
      complex(kind=8), allocatable :: bcorj(:,:) ! buffer for 
                                                 ! correlation function
c
      integer iproc ! the current processor rank
      integer nproc ! the number of processor ranks
c
      integer is, it ! counters for the correlation functions
      integer ii, ij ! counters for the natural orbitals
      integer ityp   ! type of g_ocx
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
      call ga_inquire(g_c,ityp,ii,ij)
      if (ityp.ne.MT_DCPL) then
        call errquit("wfn1_calc_ocx: wrong type of g_c",ityp,UERR)
      endif
      call ga_inquire(g_ocx,ityp,ii,ij)
      if (ityp.ne.MT_DBL) then
        call errquit("wfn1_calc_ocx: wrong type of g_ocx",ityp,UERR)
      endif
c
      call ga_distribution(g_ocx,iproc,iolo,iohi,jolo,johi)
      allocate(bocc(iolo:iohi,jolo:johi))
      allocate(boccz(iolo:iohi,jolo:johi))
      allocate(bcori(iolo:iohi,1:nmo))
      allocate(bcorj(jolo:johi,1:nmo))
      bocc = 0.0d0
      boccz = complex(0.0d0,0.0d0)
      call ga_get(g_c,iolo,iohi,1,nmo,bcori,iohi-iolo+1)
      call ga_get(g_c,jolo,johi,1,nmo,bcorj,johi-jolo+1)
      do ij = jolo, johi
        do ii = iolo, iohi
          do is = 1, ne
            do it = 1, ne
              boccz(ii,ij)
     +        = boccz(ii,ij) 
     +        - bcori(ii,icnf(is))*conjg(bcori(ii,icnf(it)))
     +        * bcorj(ij,icnf(it))*conjg(bcorj(ij,icnf(is)))
     +        * dcnt(is)*dcnt(it)
            enddo
          enddo
        enddo ! ii
      enddo ! ij
      do ij = jolo, johi
        do ii = iolo, iohi
          if (abs(aimag(boccz(ii,ij))).gt.1.0d-10) then
            write(*,*)'offending entry:',ii,ij
            do it = jolo, johi
              do is = iolo, iohi
                write(*,*)is,it,boccz(is,it)
              enddo
            enddo
            call errquit("wfn1_calc_ocx: result should be REAL "//
     +                   "instead of COMPLEX",0,UERR)
          endif
          bocc(ii,ij) = real(boccz(ii,ij))
        enddo ! ii
      enddo ! ij
      call ga_put(g_ocx,iolo,iohi,jolo,johi,bocc,iohi-iolo+1)
c
      deallocate(bcorj)
      deallocate(bcori)
      deallocate(boccz)
      deallocate(bocc)
c
      end subroutine wfn1_calc_ocx
c
c-----------------------------------------------------------------------
C>
C> \brief Add the 2-electron direct occupation numbers
C>
C> These quantities are the
C> \f$+\sum_{s,t}C_{is}C_{is}^*C_{jt}C_{jt}^*\f$ terms in the 2-electron
C> occupation numbers.
C>
C> Note that for the 2-electron density matrix the Coulomb and
C> exchange terms are separated. The reason is that the correlation
C> function Fock matrices for the Coulomb and exchange terms have to be
C> handled differently because of their nature.
C>
      subroutine wfn1_add_occd(nmo,occ,g_occ)
      implicit none
#include "global.fh"
#include "mafdecls.fh"
c
!>    The number of molecular orbitals
      integer, intent(in) :: nmo
!>    The 1-electron occupation numbers
      double precision, intent(in) :: occ(nmo)
!>    The 2-electron occupation numbers
      integer, intent(in) :: g_occ
c
      integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
      double precision, allocatable :: bocc(:,:)  ! buffer for 
                                                  ! occupation numbers
c
      integer iproc ! the current processor rank
      integer nproc ! the number of processor ranks
c
      integer ii, ij ! counters for the natural orbitals
c
      nproc = ga_nnodes()
      iproc = ga_nodeid()
c
      call ga_distribution(g_occ,iproc,iolo,iohi,jolo,johi)
      allocate(bocc(iolo:iohi,jolo:johi))
      do ij = jolo, johi
        do ii = iolo, iohi
          bocc(ii,ij) = occ(ii)*occ(ij)
        enddo
      enddo
      call ga_acc(g_occ,iolo,iohi,jolo,johi,bocc,iohi-iolo+1,1.0d0)
c
      deallocate(bocc)
c
      end subroutine wfn1_add_occd
c
c-----------------------------------------------------------------------
c
C>
C> \brief Calculate the derivative of the 2-electron exchange 
C> occupation numbers
C>
C> These quantities are the \f$-sum_s C_{js}C_{ks} terms in the
C> derivative of the 2-electron density matrix.
C>
!     subroutine wfn1_calc_docx(g_c,nmo,ne,icnf,dcnt,g_docx)
!     implicit none
!include "global.fh"
!include "mafdecls.fh"
!include "errquit.fh"
c
!     integer, intent(in) :: g_c       ! the correlation functions
!     integer, intent(in) :: nmo       ! the number of molecular orbitals
!     integer, intent(in) :: ne        ! the number of electrons
!     integer, intent(in) :: icnf(1:ne)! the electron configuration
c
!     double precision, intent(in) :: dcnt(1:ne)! the electron fraction
c
!     integer, intent(in) :: g_docx  ! the 2-electron occupation 
!                                    ! numbers derivative
c
!     integer :: iolo, iohi, jolo, johi ! block limits of occupation GA
c
!     complex(kind=8), allocatable :: bocc(:,:)  ! buffer for 
!                                                ! occupation numbers
!     complex(kind=8), allocatable :: bcori(:,:) ! buffer for 
!                                                ! correlation function
!     complex(kind=8), allocatable :: bcorj(:,:) ! buffer for 
!                                                ! correlation function
c
!     integer iproc ! the current processor rank
!     integer nproc ! the number of processor ranks
c
!     integer is     ! counters for the correlation functions
!     integer ii, ij ! counters for the natural orbitals
c
!     nproc = ga_nnodes()
!     iproc = ga_nodeid()
c
!     call ga_distribution(g_docx,iproc,iolo,iohi,jolo,johi)
!     allocate(bocc(iolo:iohi,jolo:johi))
!     allocate(bcori(iolo:iohi,1:nmo))
!     allocate(bcorj(jolo:johi,1:nmo))
!     bocc = 0.0d0
!     call ga_get(g_c,iolo,iohi,1,nmo,bcori,iohi-iolo+1)
!     call ga_get(g_c,jolo,johi,1,nmo,bcorj,johi-jolo+1)
!     do ij = jolo, johi
!       do ii = iolo, iohi
!         do is = 1, ne
!           bocc(ii,ij) = bocc(ii,ij) 
!    +                  - bcori(ii,icnf(is))*conjg(bcorj(ij,icnf(is)))
!    +                  * dcnt(is)
!         enddo
!       enddo
!     enddo
!     call ga_put(g_docx,iolo,iohi,jolo,johi,bocc,iohi-iolo+1)
c
!     deallocate(bcorj)
!     deallocate(bcori)
!     deallocate(bocc)
c
!     end subroutine wfn1_calc_docx
c
c-----------------------------------------------------------------------
