      SUBROUTINE ccsdt_1pdm_pp_mo(d_i0,d_t1,d_t2,d_t3,d_y1,d_y2,d_y3,k_i
     &0_offset,k_t1_offset,k_t2_offset,k_t3_offset,k_y1_offset,k_y2_offs
     &et,k_y3_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p1 p2 )_yt + = 1 * Sum ( h3 ) * t ( p1 h3 )_t * y ( h3 p2 )_y
C     i0 ( p1 p2 )_yt + = 1/2 * Sum ( h4 h5 p3 ) * t ( p1 p3 h4 h5 )_t * y ( h4 h5 p2 p3 )_y
C     i0 ( p1 p2 )_yt + = 1/12 * Sum ( h5 h6 h7 p3 p4 ) * t ( p1 p3 p4 h5 h6 h7 )_t * y ( h5 h6 h7 p2 p3 p4 )_y
      IMPLICIT NONE
#include "global.fh"
#include "mafdecls.fh"
#include "util.fh"
#include "errquit.fh"
#include "tce.fh"
      INTEGER d_i0
      INTEGER k_i0_offset
      INTEGER d_t1
      INTEGER k_t1_offset
      INTEGER d_y1
      INTEGER k_y1_offset
      INTEGER d_t2
      INTEGER k_t2_offset
      INTEGER d_y2
      INTEGER k_y2_offset
      INTEGER d_t3
      INTEGER k_t3_offset
      INTEGER d_y3
      INTEGER k_y3_offset
      CALL ccsdt_1pdm_pp_mo_1(d_t1,k_t1_offset,d_y1,k_y1_offset,d_i0,k_i
     &0_offset)
      CALL ccsdt_1pdm_pp_mo_2(d_t2,k_t2_offset,d_y2,k_y2_offset,d_i0,k_i
     &0_offset)
      CALL ccsdt_1pdm_pp_mo_3(d_t3,k_t3_offset,d_y3,k_y3_offset,d_i0,k_i
     &0_offset)
      RETURN
      END
      SUBROUTINE ccsdt_1pdm_pp_mo_1(d_a,k_a_offset,d_b,k_b_offset,d_c,k_
     &c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p1 p2 )_yt + = 1 * Sum ( h3 ) * t ( p1 h3 )_t * y ( h3 p2 )_y
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
      INTEGER p1b
      INTEGER p2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER h3b
      INTEGER p1b_1
      INTEGER h3b_1
      INTEGER h3b_2
      INTEGER p2b_2
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
      INTEGER l_c
      INTEGER k_c
      EXTERNAL nxtask
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p1b = noab+1,noab+nvab
      DO p2b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &).ne.4)) THEN
      IF (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+p2b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+p2b-1)) .eq. ieor(irrep_
     &y,irrep_t)) THEN
      dimc = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_1',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO h3b = 1,noab
      IF (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h3b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h3b-1)) .eq. irrep_t) TH
     &EN
      CALL TCE_RESTRICTED_2(p1b,h3b,p1b_1,h3b_1)
      CALL TCE_RESTRICTED_2(h3b,p2b,h3b_2,p2b_2)
      dim_common = int_mb(k_range+h3b-1)
      dima_sort = int_mb(k_range+p1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_1',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_1',2,MA_ERR)
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h3b_1
     & - 1 + noab * (p1b_1 - noab - 1)))
      CALL TCE_SORT_2(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+h3b-1),1,2,1.0d0)
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1',3,MA
     &_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_1',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_1',5,MA_ERR)
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p2b_2
     & - noab - 1 + nvab * (h3b_2 - 1)))
      CALL TCE_SORT_2(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h3b-1)
     &,int_mb(k_range+p2b-1),2,1,1.0d0)
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1',6,MA
     &_ERR)
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,1.0d0,dbl_mb(k_a
     &_sort),dim_common,dbl_mb(k_b_sort),dim_common,1.0d0,dbl_mb(k_c_sor
     &t),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1'
     &,7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1'
     &,8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_1',9,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+p2b-1)
     &,int_mb(k_range+p1b-1),2,1,1.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p2b -
     & noab - 1 + nvab * (p1b - noab - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1',10,M
     &A_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_1'
     &,11,MA_ERR)
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
      SUBROUTINE ccsdt_1pdm_pp_mo_2(d_a,k_a_offset,d_b,k_b_offset,d_c,k_
     &c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p1 p2 )_yt + = 1/2 * Sum ( h4 h5 p3 ) * t ( p1 p3 h4 h5 )_t * y ( h4 h5 p2 p3 )_y
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
      INTEGER p1b
      INTEGER p2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p3b
      INTEGER h4b
      INTEGER h5b
      INTEGER p1b_1
      INTEGER p3b_1
      INTEGER h4b_1
      INTEGER h5b_1
      INTEGER h4b_2
      INTEGER h5b_2
      INTEGER p2b_2
      INTEGER p3b_2
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
      INTEGER nsubh(2)
      INTEGER isubh
      INTEGER l_c
      INTEGER k_c
      DOUBLE PRECISION FACTORIAL
      EXTERNAL nxtask
      EXTERNAL FACTORIAL
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p1b = noab+1,noab+nvab
      DO p2b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &).ne.4)) THEN
      IF (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+p2b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+p2b-1)) .eq. ieor(irrep_
     &y,irrep_t)) THEN
      dimc = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_2',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p3b = noab+1,noab+nvab
      DO h4b = 1,noab
      DO h5b = h4b,noab
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p3b-1) .eq. int_mb(k_spin+h
     &4b-1)+int_mb(k_spin+h5b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+h4b-1),int_mb(k_sym+h5b-1)))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_4(p1b,p3b,h4b,h5b,p1b_1,p3b_1,h4b_1,h5b_1)
      CALL TCE_RESTRICTED_4(h4b,h5b,p2b,p3b,h4b_2,h5b_2,p2b_2,p3b_2)
      dim_common = int_mb(k_range+p3b-1) * int_mb(k_range+h4b-1) * int_m
     &b(k_range+h5b-1)
      dima_sort = int_mb(k_range+p1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_2',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_2',2,MA_ERR)
      IF ((p3b .lt. p1b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h5b_1
     & - 1 + noab * (h4b_1 - 1 + noab * (p1b_1 - noab - 1 + nvab * (p3b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p3b-1)
     &,int_mb(k_range+p1b-1),int_mb(k_range+h4b-1),int_mb(k_range+h5b-1)
     &,2,4,3,1,-1.0d0)
      END IF
      IF ((p1b .le. p3b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h5b_1
     & - 1 + noab * (h4b_1 - 1 + noab * (p3b_1 - noab - 1 + nvab * (p1b_
     &1 - noab - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p3b-1),int_mb(k_range+h4b-1),int_mb(k_range+h5b-1)
     &,1,4,3,2,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2',3,MA
     &_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_2',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_2',5,MA_ERR)
      IF ((p3b .lt. p2b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p2b_2
     & - noab - 1 + nvab * (p3b_2 - noab - 1 + nvab * (h5b_2 - 1 + noab 
     &* (h4b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h5b-1),int_mb(k_range+p3b-1),int_mb(k_range+p2b-1)
     &,4,2,1,3,-1.0d0)
      END IF
      IF ((p2b .le. p3b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p3b_2
     & - noab - 1 + nvab * (p2b_2 - noab - 1 + nvab * (h5b_2 - 1 + noab 
     &* (h4b_2 - 1)))))
      CALL TCE_SORT_4(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h4b-1)
     &,int_mb(k_range+h5b-1),int_mb(k_range+p2b-1),int_mb(k_range+p3b-1)
     &,3,2,1,4,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2',6,MA
     &_ERR)
      nsubh(1) = 1
      nsubh(2) = 1
      isubh = 1
      IF (h4b .eq. h5b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,2.0d0/FACTORIAL(
     &nsubh(1))/FACTORIAL(nsubh(2)),dbl_mb(k_a_sort),dim_common,dbl_mb(k
     &_b_sort),dim_common,1.0d0,dbl_mb(k_c_sort),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2'
     &,7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2'
     &,8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_2',9,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+p2b-1)
     &,int_mb(k_range+p1b-1),2,1,1.0d0/2.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p2b -
     & noab - 1 + nvab * (p1b - noab - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2',10,M
     &A_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_2'
     &,11,MA_ERR)
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
      SUBROUTINE ccsdt_1pdm_pp_mo_3(d_a,k_a_offset,d_b,k_b_offset,d_c,k_
     &c_offset)
C     $Id$
C     This is a Fortran77 program generated by Tensor Contraction Engine v.1.0
C     Copyright (c) Battelle & Pacific Northwest National Laboratory (2002)
C     i0 ( p1 p2 )_yt + = 1/12 * Sum ( h5 h6 h7 p3 p4 ) * t ( p1 p3 p4 h5 h6 h7 )_t * y ( h5 h6 h7 p2 p3 p4 )_y
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
      INTEGER p1b
      INTEGER p2b
      INTEGER dimc
      INTEGER l_c_sort
      INTEGER k_c_sort
      INTEGER p3b
      INTEGER p4b
      INTEGER h5b
      INTEGER h6b
      INTEGER h7b
      INTEGER p1b_1
      INTEGER p3b_1
      INTEGER p4b_1
      INTEGER h5b_1
      INTEGER h6b_1
      INTEGER h7b_1
      INTEGER h5b_2
      INTEGER h6b_2
      INTEGER h7b_2
      INTEGER p2b_2
      INTEGER p3b_2
      INTEGER p4b_2
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
      INTEGER nsuperp(2)
      INTEGER isuperp
      INTEGER nsubh(3)
      INTEGER isubh
      INTEGER l_c
      INTEGER k_c
      DOUBLE PRECISION FACTORIAL
      EXTERNAL nxtask
      EXTERNAL FACTORIAL
      nprocs = GA_NNODES()
      count = 0
      next = nxtask(nprocs,1)
      DO p1b = noab+1,noab+nvab
      DO p2b = noab+1,noab+nvab
      IF (next.eq.count) THEN
      IF ((.not.restricted).or.(int_mb(k_spin+p1b-1)+int_mb(k_spin+p2b-1
     &).ne.4)) THEN
      IF (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+p2b-1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+p2b-1)) .eq. ieor(irrep_
     &y,irrep_t)) THEN
      dimc = int_mb(k_range+p1b-1) * int_mb(k_range+p2b-1)
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c_sort,k_c_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_3',0,MA_ERR)
      CALL DFILL(dimc,0.0d0,dbl_mb(k_c_sort),1)
      DO p3b = noab+1,noab+nvab
      DO p4b = p3b,noab+nvab
      DO h5b = 1,noab
      DO h6b = h5b,noab
      DO h7b = h6b,noab
      IF (int_mb(k_spin+p1b-1)+int_mb(k_spin+p3b-1)+int_mb(k_spin+p4b-1)
     & .eq. int_mb(k_spin+h5b-1)+int_mb(k_spin+h6b-1)+int_mb(k_spin+h7b-
     &1)) THEN
      IF (ieor(int_mb(k_sym+p1b-1),ieor(int_mb(k_sym+p3b-1),ieor(int_mb(
     &k_sym+p4b-1),ieor(int_mb(k_sym+h5b-1),ieor(int_mb(k_sym+h6b-1),int
     &_mb(k_sym+h7b-1)))))) .eq. irrep_t) THEN
      CALL TCE_RESTRICTED_6(p1b,p3b,p4b,h5b,h6b,h7b,p1b_1,p3b_1,p4b_1,h5
     &b_1,h6b_1,h7b_1)
      CALL TCE_RESTRICTED_6(h5b,h6b,h7b,p2b,p3b,p4b,h5b_2,h6b_2,h7b_2,p2
     &b_2,p3b_2,p4b_2)
      dim_common = int_mb(k_range+p3b-1) * int_mb(k_range+p4b-1) * int_m
     &b(k_range+h5b-1) * int_mb(k_range+h6b-1) * int_mb(k_range+h7b-1)
      dima_sort = int_mb(k_range+p1b-1)
      dima = dim_common * dima_sort
      dimb_sort = int_mb(k_range+p2b-1)
      dimb = dim_common * dimb_sort
      IF ((dima .gt. 0) .and. (dimb .gt. 0)) THEN
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a_sort,k_a_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_3',1,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dima,'noname',l_a,k_a)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_3',2,MA_ERR)
      IF ((p4b .lt. p1b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1 + noab * (p1b_1 - noa
     &b - 1 + nvab * (p4b_1 - noab - 1 + nvab * (p3b_1 - noab - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p3b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+p1b-1),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),3,6,5,4,2,1,1.0d0)
      END IF
      IF ((p3b .lt. p1b) .and. (p1b .le. p4b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1 + noab * (p4b_1 - noa
     &b - 1 + nvab * (p1b_1 - noab - 1 + nvab * (p3b_1 - noab - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p3b-1)
     &,int_mb(k_range+p1b-1),int_mb(k_range+p4b-1),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),2,6,5,4,3,1,-1.0d0)
      END IF
      IF ((p1b .le. p3b)) THEN
      CALL GET_HASH_BLOCK(d_a,dbl_mb(k_a),dima,int_mb(k_a_offset),(h7b_1
     & - 1 + noab * (h6b_1 - 1 + noab * (h5b_1 - 1 + noab * (p4b_1 - noa
     &b - 1 + nvab * (p3b_1 - noab - 1 + nvab * (p1b_1 - noab - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_a),dbl_mb(k_a_sort),int_mb(k_range+p1b-1)
     &,int_mb(k_range+p3b-1),int_mb(k_range+p4b-1),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),1,6,5,4,3,2,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_a)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3',3,MA
     &_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b_sort,k_b_sort)) CALL
     & ERRQUIT('ccsdt_1pdm_pp_mo_3',4,MA_ERR)
      IF (.not.MA_PUSH_GET(mt_dbl,dimb,'noname',l_b,k_b)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_3',5,MA_ERR)
      IF ((p4b .lt. p2b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p2b_2
     & - noab - 1 + nvab * (p4b_2 - noab - 1 + nvab * (p3b_2 - noab - 1 
     &+ nvab * (h7b_2 - 1 + noab * (h6b_2 - 1 + noab * (h5b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),int_mb(k_range+p3b-1)
     &,int_mb(k_range+p4b-1),int_mb(k_range+p2b-1),6,3,2,1,5,4,1.0d0)
      END IF
      IF ((p3b .lt. p2b) .and. (p2b .le. p4b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p4b_2
     & - noab - 1 + nvab * (p2b_2 - noab - 1 + nvab * (p3b_2 - noab - 1 
     &+ nvab * (h7b_2 - 1 + noab * (h6b_2 - 1 + noab * (h5b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),int_mb(k_range+p3b-1)
     &,int_mb(k_range+p2b-1),int_mb(k_range+p4b-1),5,3,2,1,6,4,-1.0d0)
      END IF
      IF ((p2b .le. p3b)) THEN
      CALL GET_HASH_BLOCK(d_b,dbl_mb(k_b),dimb,int_mb(k_b_offset),(p4b_2
     & - noab - 1 + nvab * (p3b_2 - noab - 1 + nvab * (p2b_2 - noab - 1 
     &+ nvab * (h7b_2 - 1 + noab * (h6b_2 - 1 + noab * (h5b_2 - 1)))))))
      CALL TCE_SORT_6(dbl_mb(k_b),dbl_mb(k_b_sort),int_mb(k_range+h5b-1)
     &,int_mb(k_range+h6b-1),int_mb(k_range+h7b-1),int_mb(k_range+p2b-1)
     &,int_mb(k_range+p3b-1),int_mb(k_range+p4b-1),4,3,2,1,6,5,1.0d0)
      END IF
      IF (.not.MA_POP_STACK(l_b)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3',6,MA
     &_ERR)
      nsuperp(1) = 1
      nsuperp(2) = 1
      isuperp = 1
      IF (p3b .eq. p4b) THEN
      nsuperp(isuperp) = nsuperp(isuperp) + 1
      ELSE
      isuperp = isuperp + 1
      END IF
      nsubh(1) = 1
      nsubh(2) = 1
      nsubh(3) = 1
      isubh = 1
      IF (h5b .eq. h6b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      IF (h6b .eq. h7b) THEN
      nsubh(isubh) = nsubh(isubh) + 1
      ELSE
      isubh = isubh + 1
      END IF
      CALL DGEMM('T','N',dima_sort,dimb_sort,dim_common,12.0d0/FACTORIAL
     &(nsuperp(1))/FACTORIAL(nsuperp(2))/FACTORIAL(nsubh(1))/FACTORIAL(n
     &subh(2))/FACTORIAL(nsubh(3)),dbl_mb(k_a_sort),dim_common,dbl_mb(k_
     &b_sort),dim_common,1.0d0,dbl_mb(k_c_sort),dima_sort)
      IF (.not.MA_POP_STACK(l_b_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3'
     &,7,MA_ERR)
      IF (.not.MA_POP_STACK(l_a_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3'
     &,8,MA_ERR)
      END IF
      END IF
      END IF
      END DO
      END DO
      END DO
      END DO
      END DO
      IF (.not.MA_PUSH_GET(mt_dbl,dimc,'noname',l_c,k_c)) CALL ERRQUIT('
     &ccsdt_1pdm_pp_mo_3',9,MA_ERR)
      CALL TCE_SORT_2(dbl_mb(k_c_sort),dbl_mb(k_c),int_mb(k_range+p2b-1)
     &,int_mb(k_range+p1b-1),2,1,1.0d0/12.0d0)
      CALL ADD_HASH_BLOCK(d_c,dbl_mb(k_c),dimc,int_mb(k_c_offset),(p2b -
     & noab - 1 + nvab * (p1b - noab - 1)))
      IF (.not.MA_POP_STACK(l_c)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3',10,M
     &A_ERR)
      IF (.not.MA_POP_STACK(l_c_sort)) CALL ERRQUIT('ccsdt_1pdm_pp_mo_3'
     &,11,MA_ERR)
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
