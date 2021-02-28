!     submodule (wfn1) wfn1_gen_tensors
!     implicit none
!     contains
c
c-----------------------------------------------------------------------
c
!>
!> \brief calculate the 1-electron matrices
!>
!> Calculate the overlap, the 1-electron hamiltonian, and the 
!> Cauchy-Schwarz matrices.
!>
      subroutine wfn1_gen_mat(wfn1_inst,wfn1_mat)
      implicit none
!include "wfn1_param.fh"
!include "wfn1_matcnst.fh"
#include "bas.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
!>
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(inout) :: wfn1_mat
c
      integer max1e ! The integral buffer size.
      integer mem1  ! The scratch memory size for the integral code.
c
      integer l_buf, k_buf ! The buffer memory segment
      integer l_scr, k_scr ! The scratch memory segment
c
      integer iproc   ! the processor rank
      integer nproc   ! the number of processors
      integer nshell  ! the number of shells
      integer ishell  ! the i-shell number
      integer jshell  ! the j-shell number
c
      integer ilo, ihi ! the basis function limits
      integer jlo, jhi ! the basis function limits
      integer klo, khi ! the local matrix block limits
      integer llo, lhi ! the local matrix block limits
c
      character*13 pname
      pname = "wfn1_gen_mat:"
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
c
      call ga_zero(wfn1_mat%g_s)
      call ga_zero(wfn1_mat%g_tv)
      call ga_zero(wfn1_mat%g_cs)
      call ga_distribution(wfn1_mat%g_s,iproc,klo,khi,llo,lhi)
c
      if (.not.bas_numcont(wfn1_inst%basis,nshell)) call errquit
     &  (pname//" bas_numcont failed for basis",wfn1_inst%basis,
     &   BASIS_ERR)
c
      call int_mem_1e(max1e,mem1)
      if (.not. MA_push_get(MT_DBL,max1e,pname//'buf',l_buf,k_buf))
     &  call errquit(pname//' ma buf failed', max1e, MA_ERR)
      if (.not. MA_push_get(MT_DBL,mem1,pname//'scr',l_scr,k_scr))
     &  call errquit(pname//' ma scr failed', mem1, MA_ERR)
c
      do jshell = 1, nshell
        if (.not.bas_cn2bfr(wfn1_inst%basis,jshell,jlo,jhi)) then
          call errquit(pname//' bas_cn2bfr',jshell,BASIS_ERR)
        endif
        if (klo.le.jlo.and.jlo.le.khi) then
          do ishell = 1, nshell
            if (.not.bas_cn2bfr(wfn1_inst%basis,ishell,ilo,ihi)) then
              call errquit(pname//' bas_cn2bfr',ishell,BASIS_ERR)
            endif
            if (llo.le.ilo.and.ilo.le.lhi) then
c
c             Overlap
c
              call int_1eov(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_s,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Kinetic energy
c
              call int_1eke(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_tv,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Nuclear attraction integrals
c
              call int_1epe(wfn1_inst%basis,jshell,
     &                      wfn1_inst%basis,ishell,
     &                      mem1,dbl_mb(k_scr),
     &                      max1e,dbl_mb(k_buf))
              call ga_acc(wfn1_mat%g_tv,ilo,ihi,jlo,jhi,dbl_mb(k_buf),
     &                    ihi-ilo+1,1.0d0)
c
c             Cauchy-Schwarz integrals (skipped for now)
c
            endif
          enddo
        endif
      enddo
c
      if (.not. MA_pop_stack(l_scr)) call errquit
     &  (pname//' pop scr failed', 0, MA_ERR)
      if (.not. MA_pop_stack(l_buf)) call errquit
     &  (pname//' pop buf failed', 0, MA_ERR)
c
      end subroutine wfn1_gen_mat
c
c-----------------------------------------------------------------------
c
!> \brief calculate the 2-electron integral tensor
!>
!> The parallelization drives off the data distribution. If for a 
!> given processor the lowest index integral is in local memory
!> then we calculate the shell quartet.
!>
      subroutine wfn1_gen_ten2(wfn1_inst,wfn1_mat)
      implicit none
#include "bas.fh"
#include "errquit.fh"
#include "global.fh"
#include "mafdecls.fh"
!
      type(wfn1_prminst), intent(in)    :: wfn1_inst
      type(wfn1_matcnst), intent(inout) :: wfn1_mat
!
      integer :: max2e ! The integral buffer size
      integer :: mem2  ! The scratch memory size for the integral code
!
      integer :: iproc ! The processor rank
      integer :: nproc ! The number of processors
!
      integer :: mlo(4) ! lower local limits
      integer :: mhi(4) ! upper local limits
!
      integer :: blo(4) ! lower local limits
      integer :: bhi(4) ! upper local limits
      integer :: bld(4) ! leading dimensions
!
      integer :: nshell ! the number of shells in the basis set
      integer :: lsha   ! a-shell counter
      integer :: lshb   ! b-shell counter
      integer :: lshc   ! c-shell counter
      integer :: lshd   ! d-shell counter
!
      integer :: ia_lo, ia_hi ! local shell limits
      integer :: ib_lo, ib_hi ! local shell limits
      integer :: ic_lo, ic_hi ! local shell limits
      integer :: id_lo, id_hi ! local shell limits
!
      double precision, allocatable         :: scr(:)
      double precision, allocatable, target :: buf_eri(:)
      double precision, pointer             :: eri(:,:,:,:)
!
      character(len=14) pname
      pname = "wfn1_gen_ten2:"
!
      iproc = ga_nodeid()
      nproc = ga_nnodes()
!
      call int_mem_2e4c(max2e,mem2)
      allocate(scr(mem2))
      allocate(buf_eri(max2e))
!
      call ga_zero(wfn1_mat%g_2el)
      call nga_distribution(wfn1_mat%g_2el,iproc,mlo,mhi)
!
      if (.not.bas_numcont(wfn1_inst%basis,nshell)) call errquit
     &  (pname//" bas_numcont failed for basis",wfn1_inst%basis,
     &   BASIS_ERR)
!
      do lshd = 1, nshell
        if (.not. bas_cn2bfr(wfn1_inst%basis, lshd, id_lo, id_hi))
     &      call errquit('wfn1_gen_ten2: bas_cn2bfr',
     &                   wfn1_inst%basis, BASIS_ERR)
        if (id_lo.lt.mlo(4)) cycle
        if (id_lo.gt.mhi(4)) cycle
        do lshb = 1, nshell
          if (.not. bas_cn2bfr(wfn1_inst%basis, lshb, ib_lo, ib_hi))
     &        call errquit('wfn1_gen_ten2: bas_cn2bfr',
     &                     wfn1_inst%basis, BASIS_ERR)
          if (ib_lo.lt.mlo(3)) cycle
          if (ib_lo.gt.mhi(3)) cycle
          do lshc = 1, nshell
            if (.not. bas_cn2bfr(wfn1_inst%basis, lshc, ic_lo, ic_hi))
     &          call errquit('wfn1_gen_ten2: bas_cn2bfr',
     &                       wfn1_inst%basis, BASIS_ERR)
            if (ic_lo.lt.mlo(2)) cycle
            if (ic_lo.gt.mhi(2)) cycle
            do lsha = 1, nshell
              if (.not. bas_cn2bfr(wfn1_inst%basis, lsha, ia_lo, ia_hi))
     &            call errquit('wfn1_gen_ten2: bas_cn2bfr',
     &                         wfn1_inst%basis, BASIS_ERR)
              if (ia_lo.lt.mlo(1)) cycle
              if (ia_lo.gt.mhi(1)) cycle
!
!             Set integral buffer up
!
              eri(ia_lo:ia_hi,ic_lo:ic_hi,ib_lo:ib_hi,id_lo:id_hi)
     &        => buf_eri
!
!             Calculate the 2-electron integrals
!
              call int_2e4c(wfn1_inst%basis, lshd, lshb,
     &                      wfn1_inst%basis, lshc, lsha,
     &                      mem2, scr, max2e, eri)
!
!             Store the integrals
!
              blo(1) = ia_lo
              bhi(1) = ia_hi
              bld(1) = ia_hi-ia_lo+1
              blo(2) = ic_lo
              bhi(2) = ic_hi
              bld(2) = ic_hi-ic_lo+1
              blo(3) = ib_lo
              bhi(3) = ib_hi
              bld(3) = ib_hi-ib_lo+1
              blo(4) = id_lo
              bhi(4) = id_hi
              bld(4) = id_hi-id_lo+1
              call nga_put(wfn1_mat%g_2el,blo,bhi,eri,bld)
!
            enddo ! lsha
          enddo ! lshc
        enddo ! lshb
      enddo ! lshd
!
      deallocate(buf_eri)
      deallocate(scr)
!
      end subroutine wfn1_gen_ten2
c
c-----------------------------------------------------------------------
c
!     end submodule wfn1_gen_tensors
