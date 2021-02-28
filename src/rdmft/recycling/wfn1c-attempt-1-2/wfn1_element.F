c-----------------------------------------------------------------------
c
!> \brief Select the element with the largest value
!>
!> When optimizing the canonical orbitals we need to pick one matrix
!> element to polish away. Clearly we want to pick the largest one
!> to make the most progress. Subsequently we want to read the same
!> element every iteration in the line search. The later is done in 
!> the friend routine `wfn1\_fetch\_max`. Here we select the maximum
!> element and store which matrix it was found in and at what position.
!> Finally, when polishing away a redundant rotation all other rotations
!> need to be suppressed. The friend routine `wfn1\_limit\_max` 
!> changes the Fock matrices to achieve this.
!>
      function wfn1_select_max(wfn1_deriv)
      implicit none
#include "errquit.fh"
#include "wfn1_max.fh"
#include "wfn1_error.fh"
#include "wfn1_drv.fh"
!<    The result of the function
      type(wfn1_error) :: wfn1_select_max
!<    The Fock matrixes
      type(wfn1_drv), intent(in) :: wfn1_deriv
c
      type(wfn1_error) :: error
c
      character*16 :: pname
      parameter(pname="wfn1_select_max:")
c
      double precision :: val_Ca
      double precision :: val_Cb
      double precision :: val_Na
      double precision :: val_Nb
c
      integer :: pos_Ca(2)
      integer :: pos_Cb(2)
      integer :: pos_Na(2)
      integer :: pos_Nb(2)
c
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      call nga_select_elem(wfn1_deriv%g_dna,"max",val_Na,pos_Na)
      call nga_select_elem(wfn1_deriv%g_dnb,"max",val_Nb,pos_Nb)
      call nga_select_elem(wfn1_deriv%g_dca,"max",val_Ca,pos_Ca)
      call nga_select_elem(wfn1_deriv%g_dcb,"max",val_Cb,pos_Cb)
c
      if      (val_Na.eq.max(val_Na,val_Nb,val_Ca,val_Cb)) then
        mat_kind    = MAT_N
        spin_kind   = SPIN_A
        elm_pos     = pos_Na
        error%errna = val_Na
      else if (val_Nb.eq.max(val_Na,val_Nb,val_Ca,val_Cb)) then
        mat_kind    = MAT_N
        spin_kind   = SPIN_B
        elm_pos     = pos_Nb
        error%errnb = val_Nb
      else if (val_Ca.eq.max(val_Na,val_Nb,val_Ca,val_Cb)) then
        mat_kind    = MAT_C
        spin_kind   = SPIN_A
        elm_pos     = pos_Ca
        error%errca = val_Ca
      else if (val_Cb.eq.max(val_Na,val_Nb,val_Ca,val_Cb)) then
        mat_kind    = MAT_C
        spin_kind   = SPIN_B
        elm_pos     = pos_Cb
        error%errcb = val_Cb
      else
        call errquit(pname//" nothing matches",0,UERR)
      endif
c
      error%total = error%errna + error%errnb
     &            + error%errca + error%errcb
c
      wfn1_select_max = error
c
      return
      end function wfn1_select_max
c
c-----------------------------------------------------------------------
c
!> \brief Fetch the matrix element selected in `wfn1_select_max`
!>
!> During the optimization of the canonical orbitals we want to track
!> a single matrix element throughout a line search. The matrix element
!> is selected in `wfn1_select_max`. This function simply returns the
!> corresponding value.
!>
      function wfn1_fetch_max(wfn1_deriv)
      implicit none
#include "errquit.fh"
#include "wfn1_max.fh"
#include "wfn1_error.fh"
#include "wfn1_drv.fh"
      type(wfn1_error) :: wfn1_fetch_max
!<    The Fock matrices
      type(wfn1_drv), intent(in) :: wfn1_deriv
c
      type(wfn1_error) :: error
c
      character*15 :: pname
      parameter(pname="wfn1_fetch_max:")
c
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      if      (mat_kind.eq.MAT_N) then
        if      (spin_kind.eq.SPIN_A) then
          call nga_get(wfn1_deriv%g_dna,elm_pos,elm_pos,error%errna,1)
        else if (spin_kind.eq.SPIN_B) then
          call nga_get(wfn1_deriv%g_dnb,elm_pos,elm_pos,error%errnb,1)
        else
          call errquit(pname//" invalid natural spin",0,UERR)
        endif
      else if (mat_kind.eq.MAT_C) then
        if      (spin_kind.eq.SPIN_A) then
          call nga_get(wfn1_deriv%g_dca,elm_pos,elm_pos,error%errca,1)
        else if (spin_kind.eq.SPIN_B) then
          call nga_get(wfn1_deriv%g_dcb,elm_pos,elm_pos,error%errcb,1)
        else
          call errquit(pname//" invalid correlation spin",0,UERR)
        endif
      else
        call errquit(pname//" invalid matrix kind",0,UERR)
      endif
c
      error%total = error%errna + error%errnb
     &            + error%errca + error%errcb
c
      wfn1_fetch_max = error
c
      return
      end function wfn1_fetch_max
c
c-----------------------------------------------------------------------
c
!> \brief Limit the Fock matrices to the matrix element selected in
!> `wfn1_select_max`
!>
!> During the optimization of the canonical orbitals we want to track
!> a single matrix element throughout a line search. The matrix element
!> is selected in `wfn1_select_max`. When applying a redundant rotation
!> to eliminate this matrix element it is essential to suppress all
!> other rotations. This routine retrieves the relevant matrix element
!> and zeroes the rest of the Fock matrices to ensure the linesearch
!> focusses on eliminating the relevant element.
!>
      function wfn1_limit_max(wfn1_deriv)
      implicit none
#include "mafdecls.fh"
#include "errquit.fh"
#include "global.fh"
#include "wfn1_max.fh"
#include "wfn1_error.fh"
#include "wfn1_drv.fh"
      type(wfn1_error) :: wfn1_limit_max
!<    The Fock matrices
      type(wfn1_drv), intent(inout) :: wfn1_deriv
c
      type(wfn1_error) :: error
c
      character*15 :: pname
      parameter(pname="wfn1_limit_max:")
c
      integer :: elm_pos_ij(2)
      integer :: elm_pos_ji(2)
c
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
      double precision :: dmax, fmax
      integer :: ialo, iahi, jalo, jahi, nia
      integer :: iblo, ibhi, jblo, jbhi, nib
      integer :: ii, jj ! counters
c
      integer :: iproc ! this process rank
c
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      iproc = ga_nodeid()
c
      elm_pos_ij    = elm_pos
      elm_pos_ji(1) = elm_pos(2)
      elm_pos_ji(2) = elm_pos(1)
      if      (mat_kind.eq.MAT_N) then
        call ga_distribution(wfn1_deriv%g_dna,iproc,ialo,iahi,jalo,jahi)
        call ga_distribution(wfn1_deriv%g_dnb,iproc,iblo,ibhi,jblo,jbhi)
        if (ialo.ne.iblo.or.iahi.ne.ibhi.or.
     &      jalo.ne.jblo.or.jahi.ne.jbhi) then
          call errquit(pname//" mismatching dimensions natural Fock",
     &                 0,UERR)
        endif
        nia = iahi-ialo+1
        nib = ibhi-iblo+1
        allocate(bufa(ialo:iahi,jalo:jahi))
        allocate(bufb(iblo:ibhi,jblo:jbhi))
        call ga_get(wfn1_deriv%g_dna,ialo,iahi,jalo,jahi,bufa,nia)
        call ga_get(wfn1_deriv%g_dnb,iblo,ibhi,jblo,jbhi,bufb,nib)
        dmax = 0.0d0
        fmax = 0.0d0
        do jj = jalo, jahi
          do ii = ialo, iahi
            dmax = max(dmax,abs(bufa(ii,jj)-bufb(ii,jj)))
            fmax = max(fmax,abs(bufa(ii,jj)),abs(bufb(ii,jj)))
          enddo
        enddo
        deallocate(bufa)
        deallocate(bufb)
        call ga_dgop(1500,dmax,1,"max")
        call ga_dgop(1501,fmax,1,"max")
cDEBUG
c       if (iproc.eq.0) then
c         write(*,*)'HVD: wfn1_limit_max N:',dmax/fmax,dmax,fmax
c       endif
cDEBUG
        if (dmax/fmax.lt.1.0d-6.or.fmax.lt.1.0d-8) then
c
c         The alpha and beta channels are identical so this is most
c         likely a closed shell calculation and we should implement
c         the limit function to respect the equivalence.
c
          call nga_get(wfn1_deriv%g_dna,elm_pos,elm_pos,error%errna,1)
          call ga_zero(wfn1_deriv%g_dna)
          call nga_put(wfn1_deriv%g_dna,elm_pos_ij,elm_pos_ij,
     &                 error%errna,1)
          call nga_put(wfn1_deriv%g_dna,elm_pos_ji,elm_pos_ji,
     &                 -error%errna,1)
          call ga_copy(wfn1_deriv%g_dna,wfn1_deriv%g_dnb)
        else
          if      (spin_kind.eq.SPIN_A) then
            call nga_get(wfn1_deriv%g_dna,elm_pos,elm_pos,error%errna,1)
            call ga_zero(wfn1_deriv%g_dna)
            call ga_zero(wfn1_deriv%g_dnb)
            call nga_put(wfn1_deriv%g_dna,elm_pos_ij,elm_pos_ij,
     &                   error%errna,1)
            call nga_put(wfn1_deriv%g_dna,elm_pos_ji,elm_pos_ji,
     &                   -error%errna,1)
          else if (spin_kind.eq.SPIN_B) then
            call nga_get(wfn1_deriv%g_dnb,elm_pos,elm_pos,error%errnb,1)
            call ga_zero(wfn1_deriv%g_dna)
            call ga_zero(wfn1_deriv%g_dnb)
            call nga_put(wfn1_deriv%g_dnb,elm_pos_ij,elm_pos_ij,
     &                   error%errnb,1)
            call nga_put(wfn1_deriv%g_dnb,elm_pos_ji,elm_pos_ji,
     &                   -error%errnb,1)
          else
            call errquit(pname//" invalid natural spin",0,UERR)
          endif
        endif
      else if (mat_kind.eq.MAT_C) then
        call ga_distribution(wfn1_deriv%g_dca,iproc,ialo,iahi,jalo,jahi)
        call ga_distribution(wfn1_deriv%g_dcb,iproc,iblo,ibhi,jblo,jbhi)
        if (ialo.ne.iblo.or.iahi.ne.ibhi.or.
     &      jalo.ne.jblo.or.jahi.ne.jbhi) then
          call errquit(pname//" mismatching dimensions correlation Fock"
     &                 ,0,UERR)
        endif
        nia = iahi-ialo+1
        nib = ibhi-iblo+1
        allocate(bufa(ialo:iahi,jalo:jahi))
        allocate(bufb(iblo:ibhi,jblo:jbhi))
        call ga_get(wfn1_deriv%g_dca,ialo,iahi,jalo,jahi,bufa,nia)
        call ga_get(wfn1_deriv%g_dcb,iblo,ibhi,jblo,jbhi,bufb,nib)
        dmax = 0.0d0
        fmax = 0.0d0
        do jj = jalo, jahi
          do ii = ialo, iahi
            dmax = max(dmax,abs(bufa(ii,jj)-bufb(ii,jj)))
            fmax = max(fmax,abs(bufa(ii,jj)),abs(bufb(ii,jj)))
          enddo
        enddo
        deallocate(bufa)
        deallocate(bufb)
        call ga_dgop(1502,dmax,1,"max")
        call ga_dgop(1503,fmax,1,"max")
cDEBUG
c       if (iproc.eq.0) then
c         write(*,*)'HVD: wfn1_limit_max C:',dmax/fmax,dmax,fmax
c       endif
cDEBUG
        if (dmax/fmax.lt.1.0d-6.or.fmax.lt.1.0d-8) then
c
c         The alpha and beta channels are identical so this is most
c         likely a closed shell calculation and we should implement
c         the limit function to respect the equivalence.
c
          call nga_get(wfn1_deriv%g_dca,elm_pos,elm_pos,error%errca,1)
          call ga_zero(wfn1_deriv%g_dca)
          call nga_put(wfn1_deriv%g_dca,elm_pos_ij,elm_pos_ij,
     &                 error%errca,1)
          call nga_put(wfn1_deriv%g_dca,elm_pos_ji,elm_pos_ji,
     &                 -error%errca,1)
          call ga_copy(wfn1_deriv%g_dca,wfn1_deriv%g_dcb)
        else
          if      (spin_kind.eq.SPIN_A) then
            call nga_get(wfn1_deriv%g_dca,elm_pos,elm_pos,error%errca,1)
            call ga_zero(wfn1_deriv%g_dca)
            call ga_zero(wfn1_deriv%g_dcb)
            call nga_put(wfn1_deriv%g_dca,elm_pos_ij,elm_pos_ij,
     &                   error%errca,1)
            call nga_put(wfn1_deriv%g_dca,elm_pos_ji,elm_pos_ji,
     &                   -error%errca,1)
          else if (spin_kind.eq.SPIN_B) then
            call nga_get(wfn1_deriv%g_dcb,elm_pos,elm_pos,error%errcb,1)
            call ga_zero(wfn1_deriv%g_dca)
            call ga_zero(wfn1_deriv%g_dcb)
            call nga_put(wfn1_deriv%g_dcb,elm_pos_ij,elm_pos_ij,
     &                   error%errcb,1)
            call nga_put(wfn1_deriv%g_dcb,elm_pos_ji,elm_pos_ji,
     &                   -error%errcb,1)
          else
            call errquit(pname//" invalid correlation spin",0,UERR)
          endif
        endif
      else
        call errquit(pname//" invalid matrix kind",0,UERR)
      endif
c
      error%total = error%errna + error%errnb
     &            + error%errca + error%errcb
c
      wfn1_limit_max = error
c
      return
      end function wfn1_limit_max
c
c-----------------------------------------------------------------------
c
!> \brief Find and print where the Fock matrices differ most
!>
!> The parallelism seems to affect how closely the same the alpha and
!> beta Fock matrices are. This is unexpected and to trace where this
!> is coming from we need a routine to find and print where the Fock
!> matrices differ most.
!>
      subroutine wfn1_print_max(wfn1_deriv,message)
      implicit none
#include "mafdecls.fh"
#include "stdio.fh"
#include "errquit.fh"
#include "global.fh"
#include "wfn1_error.fh"
#include "wfn1_drv.fh"
!<    The Fock matrices
      type(wfn1_drv), intent(inout) :: wfn1_deriv
!<    Print message
      character*(*), intent(in) :: message
c
      character*15 :: pname
      parameter(pname="wfn1_print_max:")
c
      integer :: elm_pos_ij(2)
c
      double precision, allocatable :: bufa(:,:)
      double precision, allocatable :: bufb(:,:)
      double precision, allocatable :: dmaxg(:)
      double precision :: dmax, fmax
      integer :: ialo, iahi, jalo, jahi, nia
      integer :: iblo, ibhi, jblo, jbhi, nib
      integer :: ii, jj ! counters
c
      integer :: iroot
      integer :: iproc ! this process rank
      integer :: nproc ! the number of process ranks
c
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      call ga_sync()
c
      allocate(dmaxg(0:nproc-1))
      elm_pos_ij    = 0
      call ga_distribution(wfn1_deriv%g_dna,iproc,ialo,iahi,jalo,jahi)
      call ga_distribution(wfn1_deriv%g_dnb,iproc,iblo,ibhi,jblo,jbhi)
      if (ialo.ne.iblo.or.iahi.ne.ibhi.or.
     &    jalo.ne.jblo.or.jahi.ne.jbhi) then
        call errquit(pname//" mismatching dimensions natural Fock",
     &               0,UERR)
      endif
      nia = iahi-ialo+1
      nib = ibhi-iblo+1
      allocate(bufa(ialo:iahi,jalo:jahi))
      allocate(bufb(iblo:ibhi,jblo:jbhi))
      call ga_get(wfn1_deriv%g_dna,ialo,iahi,jalo,jahi,bufa,nia)
      call ga_get(wfn1_deriv%g_dnb,iblo,ibhi,jblo,jbhi,bufb,nib)
      dmax = 0.0d0
      fmax = 0.0d0
      do jj = jalo, jahi
        do ii = ialo, iahi
          if (abs(bufa(ii,jj)-bufb(ii,jj)).gt.dmax) then
            elm_pos_ij(1) = ii
            elm_pos_ij(2) = jj
            dmax = abs(bufa(ii,jj)-bufb(ii,jj))
            fmax = max(fmax,abs(bufa(ii,jj)),abs(bufb(ii,jj)))
          endif
        enddo
      enddo
      deallocate(bufa)
      deallocate(bufb)
      dmaxg = 0.0d0
      dmaxg(iproc) = dmax
      call ga_dgop(1504,dmaxg,nproc,"+")
      call ga_dgop(1505,dmax,1,"max")
      iroot = -1
      do ii = 0, nproc-1
        if (dmax.eq.dmaxg(ii)) then
          iroot = ii
        endif
      enddo
      call ga_brdcst(1000,fmax,ma_sizeof(MT_DBL,1,MT_BYTE),iroot)
      call ga_brdcst(1001,elm_pos_ij,ma_sizeof(MT_INT,2,MT_BYTE),
     &               iroot)
      if (iproc.eq.0) then
        write(LuOut,*)trim(message)//':N max diff: ',
     &                dmax/max(fmax,1.0d-15),fmax,elm_pos_ij
      endif
c
      elm_pos_ij    = 0
      call ga_distribution(wfn1_deriv%g_dca,iproc,ialo,iahi,jalo,jahi)
      call ga_distribution(wfn1_deriv%g_dcb,iproc,iblo,ibhi,jblo,jbhi)
      if (ialo.ne.iblo.or.iahi.ne.ibhi.or.
     &    jalo.ne.jblo.or.jahi.ne.jbhi) then
        call errquit(pname//" mismatching dimensions correlation Fock"
     &               ,0,UERR)
      endif
      nia = iahi-ialo+1
      nib = ibhi-iblo+1
      allocate(bufa(ialo:iahi,jalo:jahi))
      allocate(bufb(iblo:ibhi,jblo:jbhi))
      call ga_get(wfn1_deriv%g_dca,ialo,iahi,jalo,jahi,bufa,nia)
      call ga_get(wfn1_deriv%g_dcb,iblo,ibhi,jblo,jbhi,bufb,nib)
      dmax = 0.0d0
      fmax = 0.0d0
      do jj = jalo, jahi
        do ii = ialo, iahi
          if (abs(bufa(ii,jj)-bufb(ii,jj)).gt.dmax) then
            elm_pos_ij(1) = ii
            elm_pos_ij(2) = jj
            dmax = abs(bufa(ii,jj)-bufb(ii,jj))
            fmax = max(fmax,abs(bufa(ii,jj)),abs(bufb(ii,jj)))
          endif
        enddo
      enddo
      deallocate(bufa)
      deallocate(bufb)
      dmaxg = 0.0d0
      dmaxg(iproc) = dmax
      call ga_dgop(1506,dmaxg,nproc,"+")
      call ga_dgop(1507,dmax,1,"max")
      iroot = -1
      do ii = 0, nproc-1
        if (dmax.eq.dmaxg(ii)) then
          iroot = ii
        endif
      enddo
      call ga_brdcst(1002,fmax,ma_sizeof(MT_DBL,1,MT_BYTE),iroot)
      call ga_brdcst(1003,elm_pos_ij,ma_sizeof(MT_INT,2,MT_BYTE),
     &               iroot)
      if (iproc.eq.0) then
        write(LuOut,*)trim(message)//':C max diff: ',
     &                dmax/max(fmax,1.0d-15),fmax,elm_pos_ij
      endif
      deallocate(dmaxg)
c
      return
      end subroutine wfn1_print_max
c
c-----------------------------------------------------------------------
