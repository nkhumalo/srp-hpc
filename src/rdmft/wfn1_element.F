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
