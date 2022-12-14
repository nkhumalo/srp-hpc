      SUBROUTINE ccsdtq_lr_alpha_15_50_2(d_a,k_a_offset,d_b,k_b_offset,d
     &_c,k_c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i2 ( p13 h16 )_ytrbt + = 1/12 * Sum ( h9 h10 h11 h12 p5 p6 p7 ) * t ( p5 p6 p7 p13 h9 h10 h11 h12 )_t * i3 ( h9 h10 h11 h12 h16 p5 p6 p7 )_ytrb
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_a
      INTEGER k_a_offset
      INTEGER d_b
      INTEGER k_b_offset
      INTEGER d_c
      INTEGER k_c_offset
      INTEGER nxtask
      INTEGER next
      INTEGER nprocs
      INTEGER count
      INTEGER p13b
      INTEGER h16b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p5b
      INTEGER p6b
      INTEGER p7b
      INTEGER h9b
      INTEGER h10b
      INTEGER h11b
      INTEGER h12b
      INTEGER p13b_1
      INTEGER p5b_1
      INTEGER p6b_1
      INTEGER p7b_1
      INTEGER h9b_1
      INTEGER h10b_1
      INTEGER h11b_1
      INTEGER h12b_1
      INTEGER h9b_2
      INTEGER h10b_2
      INTEGER h11b_2
      INTEGER h12b_2
      INTEGER h16b_2
      INTEGER p5b_2
      INTEGER p6b_2
      INTEGER p7b_2
      INTEGER dim_common
      INTEGER dima_sort
      INTEGER dima
      INTEGER dimb_sort
      INTEGER dimb
      INTEGER l_a_sort
      INTEGER k_a_sort
      INTEGER l_a
      INTEGER k_a
      INTEGER l_b_sort
      INTEGER k_b_sort
      INTEGER l_b
      INTEGER k_b
      INTEGER nsuperp(3)
      INTEGER isuperp
      INTEGER nsubh(4)
      INTEGER isubh
      INTEGER l_c
      INTEGER k_c
      DOUBLE PRECISION FACTORIAL
      EXTERNAL nxtask
      EXTERNAL FACTORIAL
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p13b = noab+1,noab+nvab
      DO h16b = 1,noab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p13b-1)+int_mb(k_spin+h16b
     &-1).ne.4)) THEN
      IF (int_mb(k_spin+p13b-1) .eq. int_mb(k_spin+h16b-1)) THEN
      IF (ieor(int_mb(k_sym+p13b-1),int_mb(k_sym+h16b-1)) .eq. ieor(irre
     &p_y,ieor(irrep_trb,irrep_t))) THEN
      dimc = int_mb(k_range+p13b-1) * int_mb(k_range+h16b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdtq_lr_alpha_15_50_2',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p5b = noab+1,noab+nvab
      DO p6b = p5b,noab+nvab
      DO p7b = p6b,noab+nvab
      DO h9b = 1,noab
      DO h10b = h9b,noab
      DO h11b = h10b,noab
      DO h12b = h11b,noab
      IF (int_mb(k_spin+p13b-1)+int_mb(k_spin+p5b-1)+int_mb(k_spin+p6b-1
     &)+int_mb(k_spin+p7b-1) .eq. int_mb(k_spin+h9b-1)+int_mb(k_spin+h10
     &b-1)+int_mb(k_spin+h11b-1)+int_mb(k_spin+h12b-1)) THEN
      IF (ieor(int_mb(k_sym+p13b-1),ieor(int_mb(k_sym+p5b-1),ieor(int_mb
     &(k_sym+p6b-1),ieor(int_mb(k_sym+p7b-1),ieor(int_mb(k_sym+h9b-1),ie
     &or(int_mb(k_sym+h10b-1),ieor(int_mb(k_sym+h11b-1),int_mb(k_sym+h12
     &b-1)))))))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_8(p13b,p5b,p6b,p7b,h9b,h10b,h11b,h12b,p13b_1,p
     &5b_1,p6b_1,p7b_1,h9b_1,h10b_1,h11b_1,h12b_1)
      CALL TCE_RESTRICTED_8(h9b,h10b,h11b,h12b,h16b,p5b,p6b,p7b,h9b_2,h1
     &0b_2,h11b_2,h12b_2,h16b_2,p5b_2,p6b_2,p7b_2)
      dim_common = int_mb(k_range+p5b-1) * int_mb(k_range+p6b-1) * int_m
     &b(k_range+p7b-1) * int_mb(k_range+h9b-1) * int_mb(k_range+h10b-1) 
     &* int_mb(k_range+h11b-1) * int_mb(k_range+h12b-1)
      dima_sort = int_mb(k_range+p13b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+h16b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdtq_lr_alpha_15_50_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdtq_lr_alpha_15_50_2',2,MA_ERR)
      IF ((p7b .le. p13b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h12b_
     &1 - 1 + noab * (h11b_1 - 1 + noab * (h10b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (p13b_1 - noab - 1 + nvab * (p7b_1 - noab - 1 + nvab * 
     &(p6b_1 - noab - 1 + nvab * (p5b_1 - noab - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+p7b-1),int_mb(k_range+p13b-1
     &),int_mb(k_range+h9b-1),int_mb(k_range+h10b-1),int_mb(k_range+h11b
     &-1),int_mb(k_range+h12b-1),4,8,7,6,5,3,2,1,1.0d0)
      END IF
      IF ((p6b .le. p13b) .and. (p13b .lt. p7b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h12b_
     &1 - 1 + noab * (h11b_1 - 1 + noab * (h10b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (p7b_1 - noab - 1 + nvab * (p13b_1 - noab - 1 + nvab * 
     &(p6b_1 - noab - 1 + nvab * (p5b_1 - noab - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p6b-1),int_mb(k_range+p13b-1),int_mb(k_range+p7b-1
     &),int_mb(k_range+h9b-1),int_mb(k_range+h10b-1),int_mb(k_range+h11b
     &-1),int_mb(k_range+h12b-1),3,8,7,6,5,4,2,1,-1.0d0)
      END IF
      IF ((p5b .le. p13b) .and. (p13b .lt. p6b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h12b_
     &1 - 1 + noab * (h11b_1 - 1 + noab * (h10b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (p7b_1 - noab - 1 + nvab * (p6b_1 - noab - 1 + nvab * (
     &p13b_1 - noab - 1 + nvab * (p5b_1 - noab - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p5b-1)
     &,int_mb(k_range+p13b-1),int_mb(k_range+p6b-1),int_mb(k_range+p7b-1
     &),int_mb(k_range+h9b-1),int_mb(k_range+h10b-1),int_mb(k_range+h11b
     &-1),int_mb(k_range+h12b-1),2,8,7,6,5,4,3,1,1.0d0)
      END IF
      IF ((p13b .lt. p5b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h12b_
     &1 - 1 + noab * (h11b_1 - 1 + noab * (h10b_1 - 1 + noab * (h9b_1 - 
     &1 + noab * (p7b_1 - noab - 1 + nvab * (p6b_1 - noab - 1 + nvab * (
     &p5b_1 - noab - 1 + nvab * (p13b_1 - noab - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p13b-1
     &),int_mb(k_range+p5b-1),int_mb(k_range+p6b-1),int_mb(k_range+p7b-1
     &),int_mb(k_range+h9b-1),int_mb(k_range+h10b-1),int_mb(k_range+h11b
     &-1),int_mb(k_range+h12b-1),1,8,7,6,5,4,3,2,-1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdtq_lr_alpha_15_50_2'
     &,3,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdtq_lr_alpha_15_50_2',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdtq_lr_alpha_15_50_2',5,MA_ERR)
      IF ((h16b .le. p5b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p7b_2
     & - noab - 1 + nvab * (p6b_2 - noab - 1 + nvab * (p5b_2 - noab - 1 
     &+ nvab * (h16b_2 - 1 + noab * (h12b_2 - 1 + noab * (h11b_2 - 1 + n
     &oab * (h10b_2 - 1 + noab * (h9b_2 - 1)))))))))
      CALL TCE_SORT_8(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h9b-1)
     &,int_mb(k_range+h10b-1),int_mb(k_range+h11b-1),int_mb(k_range+h12b
     &-1),int_mb(k_range+h16b-1),int_mb(k_range+p5b-1),int_mb(k_range+p6
     &b-1),int_mb(k_range+p7b-1),5,4,3,2,1,8,7,6,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdtq_lr_alpha_15_50_2'
     &,6,MA_ERR)
      nsuperp(1) = 1
      nsuperp(2) = 1
      nsuperp(3) = 1
      isuperp = 1
      IF (p5b .eq. p6b) THEN
      nsuperp(isuperp) = nsuperp(isuperp) + 1
      ELSE
      isuperp = isuperp + 1
      END IF
      IF (p6b .eq. p7b) THEN
      nsuperp(isuperp) = nsuperp(isuperp) + 1
      ELSE
      isuperp = isuperp + 1
      END IF
      nsubh(1) = 1
      nsubh(2) = 1
      nsubh(3) = 1
      nsubh(4) = 1
      isubh = 1
      IF (h9b .eq. h10b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      IF (h10b .eq. h11b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      IF (h11b .eq. h12b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,144.0d0/FACTORIA
     &L(nsuperp(1))/FACTORIAL(nsuperp(2))/FACTORIAL(nsuperp(3))/FACTORIA
     &L(nsubh(1))/FACTORIAL(nsubh(2))/FACTORIAL(nsubh(3))/FACTORIAL(nsub
     &h(4)),dbl_mb(k_a_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d
     &0,dbl_mb(k_c_sort),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdtq_lr_alpha_15_
     &50_2',7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdtq_lr_alpha_15_
     &50_2',8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdtq_lr_alpha_15_50_2',9,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+h16b-1
     &),int_mb(k_range+p13b-1),2,1,1.0d0/12.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(h16b 
     &- 1 + noab * (p13b - noab - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdtq_lr_alpha_15_50_2'
     &,10,MA_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdtq_lr_alpha_15_
     &50_2',11,MA_ERR)
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
