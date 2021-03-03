c-----------------------------------------------------------------------
c
!> \brief Evaluate the 2-electron energy based on the wavefunction
!>        model
!>
      subroutine wfn1_e_g_2el_wfn(wfn1_param,wfn1_mat,wfn1_wave,
     &                            g_bfbf,g_4mo,g_4t1,E_2el)
      implicit none
#include "errquit.fh"
#include "mafdecls.fh"
#include "global.fh"
!>    Constant parameters for this calculation
      type(wfn1_prm),     intent(in)    :: wfn1_param
!>    Matrix constants for this calculation
      type(wfn1_matcnst), intent(in)    :: wfn1_mat
!>    The wave function
      type(wfn1_wfn),     intent(in)    :: wfn1_wave
!>    The energy terms
      double precision,   intent(inout) :: E_2el
!>    GA handle for a 2-dimensional temporary array
      integer,            intent(in)    :: g_bfbf
!>    GA handle for the 2-electron integrals in MO basis
      integer,            intent(in)    :: g_4mo
!>    GA handle for a 4-dimensional temporary array
      integer,            intent(in)    :: g_4t1
c
c     Local variables
c
      integer :: nea   ! number of alpha-electrons
      integer :: neb   ! number of beta-electrons
      integer :: nmo   ! number of molecular orbitals
      integer :: iproc ! this processor rank
      integer :: nproc ! the number of ranks
c
      integer :: ilo(4), ihi(4), ild(4)
      integer :: il1, ih1, ld1
      integer :: il2, ih2, ld2
      integer :: il3, ih3, ld3
      integer :: il4, ih4, ld4
c
      integer :: ii, ij, ik, il
      integer :: ir, is
c
      integer :: indxl, indxh
c
      complex(kind=8), allocatable :: cai(:,:)
      complex(kind=8), allocatable :: caj(:,:)
      complex(kind=8), allocatable :: cak(:,:)
      complex(kind=8), allocatable :: cal(:,:)
      complex(kind=8), allocatable :: cbi(:,:)
      complex(kind=8), allocatable :: cbj(:,:)
      complex(kind=8), allocatable :: cbk(:,:)
      complex(kind=8), allocatable :: cbl(:,:)
      complex(kind=8), allocatable :: daij(:,:)
      complex(kind=8), allocatable :: dbkl(:,:)
c
      complex(kind=8), parameter :: zero = complex(0.0d0,0.0d0)
c
      complex(kind=8)           :: t_2el
      double precision, pointer :: eri(:,:,:,:)
      target :: dbl_mb
c
c     Code
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      nea   = wfn1_param%nea
      neb   = wfn1_param%neb
      nmo   = wfn1_param%nmo
      t_2el = zero
c
c     Do alpha-beta electron interaction
c
      call wfn1_4index_ab(wfn1_param,wfn1_wave%g_na,wfn1_wave%g_nb,
     &                    wfn1_mat%g_2el,g_4mo,g_4t1,g_bfbf)
cDEBUG
c     call util_flush(6)
c     call ga_print(wfn1_mat%g_2el)
c     call ga_print(g_4mo)
c     call util_flush(6)
c     call errquit("stop here")
cDEBUG
      call nga_distribution(g_4mo,iproc,ilo,ihi)
      il1 = ilo(1)
      il2 = ilo(2)
      il3 = ilo(3)
      il4 = ilo(4)
      ih1 = ihi(1)
      ih2 = ihi(2)
      ih3 = ihi(3)
      ih4 = ihi(4)
      ld1 = ih1-il1+1
      ld2 = ih2-il2+1
      ld3 = ih3-il3+1
      ld4 = ih4-il4+1
      ild(1) = ld1
      ild(2) = ld2
      ild(3) = ld3
      ild(4) = ld4
      allocate(cai(il1:ih1,1:nea),caj(il2:ih2,1:nea))
      allocate(cbk(il3:ih3,1:neb),cbl(il4:ih4,1:neb))
      allocate(daij(il1:ih1,il2:ih2),dbkl(il3:ih3,il4:ih4))
      call ga_get(wfn1_wave%g_ca,il1,ih1,1,nea,cai,ld1)
      call ga_get(wfn1_wave%g_ca,il2,ih2,1,nea,caj,ld2)
      call ga_get(wfn1_wave%g_cb,il3,ih3,1,neb,cbk,ld3)
      call ga_get(wfn1_wave%g_cb,il4,ih4,1,neb,cbl,ld4)
      daij = zero
      dbkl = zero
      do ir = 1, nea
        do ij = il2, ih2
          do ii = il1, ih1
            daij(ii,ij) = daij(ii,ij) + cai(ii,ir)*conjg(caj(ij,ir))
          enddo ! ii
        enddo ! ij
      enddo ! ir
      do ir = 1, neb
        do il = il4, ih4
          do ik = il3, ih3
            dbkl(ik,il) = dbkl(ik,il) + cbk(ik,ir)*conjg(cbl(il,ir))
          enddo ! il
        enddo ! ik
      enddo ! ir
      call nga_access(g_4mo,ilo,ihi,indxl,ild)
      indxh = indxl + ild(1)*ild(2)*ild(3)*ild(4)-1
      eri(il1:ih1,il2:ih2,il3:ih3,il4:ih4) => dbl_mb(indxl:indxh)
      do il = il4, ih4
        do ik = il3, ih3
          do ij = il2, ih2
            do ii = il1, ih1
              if (ii.eq.ik.and.ij.eq.il) then
                t_2el = t_2el + eri(ii,ij,ik,il)*daij(ii,ij)*dbkl(ik,il)
cDEBUG
!     write(*,'(4i4,f12.6,3(f10.6,SP,f8.6,"i "))')ii,ij,ik,il,
!    &eri(ii,ij,ik,il),daij(ii,ij),dbkl(ik,il),daij(ii,ij)*dbkl(ik,il)
cDEBUG
              endif
            enddo ! ii
          enddo ! ij
        enddo ! ik
      enddo ! il
      deallocate(daij,dbkl)
      deallocate(cbk,cbl)
      deallocate(cai,caj)
c
c     Do alpha-alpha electron interaction
c
      call wfn1_4index_aa(wfn1_param,wfn1_wave%g_na,
     &                    wfn1_mat%g_2el,g_4mo,g_4t1,g_bfbf)
c
c     Do beta-beta electron interaction
c
      call wfn1_4index_aa(wfn1_param,wfn1_wave%g_nb,
     &                    wfn1_mat%g_2el,g_4mo,g_4t1,g_bfbf)
c
      call ga_zgop(1302,t_2el,1,"+")
      if (abs(aimag(t_2el)).gt.1.0d-10) then
        write(*,*)"Energy = ",t_2el
        call errquit("wfn1_e_g_2el_wfn: energy not real",0,UERR)
      endif
cDEBUG
c       call errquit("wfn1_e_g_2el_wfn: energy not real",0,UERR)
cDEBUG
      E_2el = 2.0d0*real(t_2el)
c
      end subroutine wfn1_e_g_2el_wfn
c
c-----------------------------------------------------------------------
