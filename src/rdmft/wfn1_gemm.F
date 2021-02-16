c-----------------------------------------------------------------------
c
!> \brief A matrix-matrix multiplication routine for real or complex
!> matrices
!>
!> We need to be able to do matrix-matrix multiplications on either
!> REAL or COMPLEX matrices. This routine calls the appropriate
!> global array routines based on the type of the matrices.
!> The parameters alpha and beta are typically simple scale factors, so
!> we can accept them as REAL and convert them as appropriate.
!>
      subroutine wfn1_gemm(ta,tb,m,n,k,alpha,g_a,g_b,beta,g_c)
      implicit none
#include "mafdecls.fh"
#include "global.fh"
#include "errquit.fh"
      character(len=1), intent(in) :: ta
      character(len=1), intent(in) :: tb
      integer, intent(in) :: m
      integer, intent(in) :: n
      integer, intent(in) :: k
      real(kind=8), intent(in) :: alpha
      real(kind=8), intent(in) :: beta
      integer, intent(in) :: g_a
      integer, intent(in) :: g_b
      integer, intent(in) :: g_c
c
      integer :: itype_a, itype_b, itype_c
      integer :: ii, jj
c
      call ga_inquire(g_a,itype_a,ii,jj)
      call ga_inquire(g_b,itype_b,ii,jj)
      call ga_inquire(g_c,itype_c,ii,jj)
c
      if (itype_a.ne.itype_b.or.itype_a.ne.itype_c) then
        if (itype_a.eq.MT_DBL) then
          write(*,*)"g_a is DOUBLE"
        else if (itype_a.eq.MT_DCPL) then
          write(*,*)"g_a is COMPLEX"
        else
          write(*,*)"g_a is ",itype_a
        endif
        if (itype_b.eq.MT_DBL) then
          write(*,*)"g_b is DOUBLE"
        else if (itype_b.eq.MT_DCPL) then
          write(*,*)"g_b is COMPLEX"
        else
          write(*,*)"g_b is ",itype_b
        endif
        if (itype_c.eq.MT_DBL) then
          write(*,*)"g_c is DOUBLE"
        else if (itype_c.eq.MT_DCPL) then
          write(*,*)"g_c is COMPLEX"
        else
          write(*,*)"g_c is ",itype_c
        endif
        call errquit("wfn1_gemm: not all matrix types are the same",0,
     &               MA_ERR)
      endif
c
      if (itype_a.eq.MT_DBL) then
        call ga_dgemm(ta,tb,m,n,k,alpha,g_a,g_b,beta,g_c)
      else if (itype_a.eq.MT_DCPL) then
        call ga_zgemm(ta,tb,m,n,k,complex(alpha,0.0d0),g_a,g_b,
     &                            complex(beta,0.0d0),g_c)
      else
        call errquit("wfn1_gemm: invalid type",itype_a,MA_ERR)
      endif
      end subroutine wfn1_gemm
c
c-----------------------------------------------------------------------
