      subroutine tce_mrcc_jacobi_t1(d_r1,d_t1,k_t1_offset,iter,iref)

      implicit none
#include "global.fh"
#include "mafdecls.fh"
#include "sym.fh"
#include "util.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "tce.fh"
#include "tce_main.fh"
#include "tce_mrcc.fh"

      integer d_r1
      integer d_t1
      integer p1b
      integer h2b
      integer p1
      integer h2
      integer k_t1_offset
      integer size
      integer l_r1,k_r1
      integer l_t1,k_t1
      integer i
      integer nprocs
      integer count
      integer next
cc      integer nxtask
c
      integer iter
      integer iref
      double precision denom,cit
      double precision d1,d2
      integer orbspin1,orbspin2
      integer orbindex1,orbindex2
      integer blck1,blck2
      integer off1,off2
c
cc      external nxtask
      INTEGER NXTASK
      EXTERNAL NXTASK
      logical nodezero
      logical noloadbalance
c *** shift ***
      double precision shift
c *************
c
      if(iter.le.100) then
       shift=-zlshift
      else
       shift=0.0d0
      end if
c
c     ================
c     Loop over blocks
c     ================
c
      nodezero = (ga_nodeid().eq.0)
      noloadbalance = ((ioalg.eq.4).or.
     1                ((ioalg.eq.6).and.(.not.fileisga(d_r1))))
      nprocs = ga_nnodes()
      count = 0
cc      next = nxtask(nprocs,1)
      call ga_sync()
      next = NXTASK(nprocs, 1)
      do p1b = noab+1,noab+nvab
        do h2b = 1,noab
          if (noloadbalance.or.(next.eq.count)) then
            if (int_mb(k_spin+p1b-1) .eq. int_mb(k_spin+h2b-1)) then 
            if ((.not.restricted).or.(int_mb(k_spin+p1b-1)
     1        +int_mb(k_spin+h2b-1).ne.4)) then
            if (ieor(int_mb(k_sym+p1b-1),int_mb(k_sym+h2b-1)).eq.0) then
              size = int_mb(k_range+p1b-1) * int_mb(k_range+h2b-1)
              if (.not.ma_push_get(mt_dbl,size,'r1',l_r1,k_r1))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)
              call get_hash_block(d_r1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
c           if(lsubterm) then
              if (.not.ma_push_get(mt_dbl,size,'t1',l_t1,k_t1))
     1          call errquit('tce_jacobi_t1: MA problem',0,MA_ERR)
              call get_hash_block(d_t1,dbl_mb(k_t1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
c           endif
              i = 0
              do p1 = 1,int_mb(k_range+p1b-1)
                do h2 = 1,int_mb(k_range+h2b-1)
                  i = i + 1

       denom = (-dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
     2              +dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1)
     3              +mrccshift)

c       if((abs(denom).lt.0.01).and.
c     1 (abs(dbl_mb(k_r1+i-1)/denom).gt.1.0d0))
c     1  then
c       write(6,"('1DENOM CLOSE TO ZERO: ',4F16.12,4F16.12)")denom,
c     1 dbl_mb(k_r1+i-1)
c       endif

       if(lsubterm) then
        cit = dbl_mb(k_r1+i-1)-(mrccshift*dbl_mb(k_t1+i-1))
       else
        cit = dbl_mb(k_r1+i-1)
       endif

       if(lusesamefock_it) then
         orbspin1 = int_mb(k_spin+p1b-1)-1
         orbspin2 = int_mb(k_spin+h2b-1)-1

         orbindex1 = (1-orbspin1+int_mb(k_mo_indexm(iref)+
     1 int_mb(k_offset+p1b-1)+p1-1))/2
         orbindex2 = (1-orbspin2+int_mb(k_mo_indexm(iref)+
     1 int_mb(k_offset+h2b-1)+h2-1))/2

       blck1 = orbinblck(orbindex1,orbspin1+1,1)
       blck2 = orbinblck(orbindex2,orbspin2+1,1)

       off1 = offsetinblck(orbindex1,orbspin1+1,1)
       off2 = offsetinblck(orbindex2,orbspin2+1,1)

      d1 = -dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)+blck1-1)+off1)
      d2 = dbl_mb(k_evl_sortedm(1)+int_mb(k_offsetm(1)+blck2-1)+off2)

c       write(6,"('T1: ',4F16.12)")d1,
c     1 -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1),
c     2 d2,
c     2 dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1)

       else
        d1 = -dbl_mb(k_evl_sorted+int_mb(k_offset+p1b-1)+p1-1)
        d2 = dbl_mb(k_evl_sorted+int_mb(k_offset+h2b-1)+h2-1)
       endif

       dbl_mb(k_r1+i-1) = cit/(d1+d2+shift+mrccshift)

c       if(abs(dbl_mb(k_r1+i-1)).gt.20.0d0) then
c           if(dbl_mb(k_r1+i-1).lt.0.0d0) then
c               dbl_mb(k_r1+i-1)=-20.0d0
c           else
c               dbl_mb(k_r1+i-1)=20.0d0
c           endif
c           write(6,"('1RESIDUE HAS BEEN CUTED')")
c       endif
      if(iter .lt. 4) then
       if(.not. lreadt) then
        if(abs(dbl_mb(k_r1+i-1)).gt. 0.1d0) then
         if(dbl_mb(k_r1+i-1).lt.0.0d0) then
               dbl_mb(k_r1+i-1)=-0.01d0
         else
               dbl_mb(k_r1+i-1)=0.01d0
         endif
c        if(nodezero) 
c          write(luout,*)"t1 residue has been modified to 0.01"
c        if(nodezero) call util_flush(LuOut)
        endif
       endif
      endif

                enddo
              enddo
              call add_hash_block(d_t1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
c update of the res.-single vector to the form of increment used in DIIS proc.
c              call put_hash_block(d_r1,dbl_mb(k_r1),size,
c     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
cc
              if (nodezero.and.util_print('t1',print_debug)) then
                call get_hash_block(d_t1,dbl_mb(k_r1),size,
     1          int_mb(k_t1_offset),((p1b-noab-1)*noab+h2b-1))
                call ma_print_compact(dbl_mb(k_r1),size,1,'t1')
              endif
c           if(lsubterm) then
              if (.not.ma_pop_stack(l_t1))
     1          call errquit('tce_jacobi_t1: MA problem',1,MA_ERR)
c           endif
              if (.not.ma_pop_stack(l_r1))
     1          call errquit('tce_jacobi_t1: MA problem',1,MA_ERR)
            endif
            endif
            endif
cc            next = nxtask(nprocs,1)
      next = NXTASK(nprocs, 1)
          endif
          count = count + 1
        enddo
      enddo
cc      next = nxtask(-nprocs,1)
      next = NXTASK(-nprocs, 1)
      call ga_sync()
      call util_flush(LuOut)
      return
      end
c
c
c

c $Id$
