      SUBROUTINE cr_eomccsd_t_N2_1_7_1(d_a,k_a_offset,d_c,k_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( h12 p7 )_f + = 1 * f ( h12 p7 )_f
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER h12b
      INTEGER p7b
      INTEGER dimc
      INTEGER h12b_1
      INTEGER p7b_1
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO h12b = 1,noab
      DO p7b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+h12b-1)+int_mb(k_spin+p7b-
     &1).ne.4)) THEN
      IF (int_mb(k_spin+h12b-1) .eq. int_mb(k_spin+p7b-1)) THEN
      IF (ieor(int_mb(k_sym+h12b-1),int_mb(k_sym+p7b-1)) .eq. irrep_f) T
     &HEN
      dimc = int_mb(k_range+h12b-1) * int_mb(k_range+p7b-1)
      CALL TCE_RESTRICTED_2(h12b,p7b,h12b_1,p7b_1)
      dim_common = 1
      dima_sort = int_mb(k_range+h12b-1) * int_mb(k_range+p7b-1)
      dima = dim_common * dima_sort
      IF (dima .gt. 0) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('cr_eomccsd_t_N2_1_7_1',0,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &cr_eomccsd_t_N2_1_7_1',1,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(p7b_1
     & - 1 + (noab+nvab) * (h12b_1 - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+h12b-1
     &),int_mb(k_range+p7b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('cr_eomccsd_t_N2_1_7_1',2
     &,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &cr_eomccsd_t_N2_1_7_1',3,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_a_sort),dbl_mb(k_c),int_mb(k_range+p7b-1)
     &,int_mb(k_range+h12b-1),2,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p7b -
     & noab - 1 + nvab * (h12b - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('cr_eomccsd_t_N2_1_7_1',4
     &,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('cr_eomccsd_t_N2_1_7
     &_1',5,MA_ERR)
      END IF
      END IF
      END IF
      END IF
      next = nxtask(nprocs,1)
      END IF
      count = count + 1
      END DO
      END DO
      next = nxtask(-nprocs,1)
      call GA_SYNC()
      RETURN
      END
