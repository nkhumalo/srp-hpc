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
#include "errquit.fh"
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
      error%total = 0.0d0
      error%errna = 0.0d0
      error%errnb = 0.0d0
      error%errca = 0.0d0
      error%errcb = 0.0d0
c
      elm_pos_ij    = elm_pos
      elm_pos_ji(1) = elm_pos(2)
      elm_pos_ji(2) = elm_pos(1)
      if      (mat_kind.eq.MAT_N) then
        if      (spin_kind.eq.SPIN_A) then
          call nga_get(wfn1_deriv%g_dna,elm_pos,elm_pos,error%errna,1)
          call ga_zero(wfn1_deriv%g_dna)
          call ga_zero(wfn1_deriv%g_dnb)
          call nga_put(wfn1_deriv%g_dna,elm_pos_ij,elm_pos_ij,
     &                 error%errna,1)
          call nga_put(wfn1_deriv%g_dna,elm_pos_ji,elm_pos_ji,
     &                 -error%errna,1)
        else if (spin_kind.eq.SPIN_B) then
          call nga_get(wfn1_deriv%g_dnb,elm_pos,elm_pos,error%errnb,1)
          call ga_zero(wfn1_deriv%g_dna)
          call ga_zero(wfn1_deriv%g_dnb)
          call nga_put(wfn1_deriv%g_dnb,elm_pos_ij,elm_pos_ij,
     &                 error%errnb,1)
          call nga_put(wfn1_deriv%g_dnb,elm_pos_ji,elm_pos_ji,
     &                 -error%errnb,1)
        else
          call errquit(pname//" invalid natural spin",0,UERR)
        endif
      else if (mat_kind.eq.MAT_C) then
        if      (spin_kind.eq.SPIN_A) then
          call nga_get(wfn1_deriv%g_dca,elm_pos,elm_pos,error%errca,1)
          call ga_zero(wfn1_deriv%g_dca)
          call ga_zero(wfn1_deriv%g_dcb)
          call nga_put(wfn1_deriv%g_dca,elm_pos_ij,elm_pos_ij,
     &                 error%errca,1)
          call nga_put(wfn1_deriv%g_dca,elm_pos_ji,elm_pos_ji,
     &                 -error%errca,1)
        else if (spin_kind.eq.SPIN_B) then
          call nga_get(wfn1_deriv%g_dcb,elm_pos,elm_pos,error%errcb,1)
          call ga_zero(wfn1_deriv%g_dca)
          call ga_zero(wfn1_deriv%g_dcb)
          call nga_put(wfn1_deriv%g_dcb,elm_pos_ij,elm_pos_ij,
     &                 error%errcb,1)
          call nga_put(wfn1_deriv%g_dcb,elm_pos_ji,elm_pos_ji,
     &                 -error%errcb,1)
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
      wfn1_limit_max = error
c
      return
      end function wfn1_limit_max
c
c-----------------------------------------------------------------------
