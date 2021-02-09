c-----------------------------------------------------------------------
      module wfn1
      implicit none
#include "wfn1_param.fh"
#include "wfn1_wfn.fh"
#include "wfn1_error.fh"
#include "wfn1_matcnst.fh"
#include "wfn1_en.fh"
#include "wfn1_drv.fh"
      interface 
        module subroutine wfn1_input(rtdb)
          implicit none
          integer, intent (inout) :: rtdb
        end subroutine
        module function wfn1_task_energy(rtdb) result(x)
          implicit none
          integer, intent(inout) :: rtdb
          logical :: x
        end function
        module subroutine wfn1_allocate_mat(wfn1_param,wfn1_inst,
     &                                      wfn1_mat)
          implicit none
          type(wfn1_prm),     intent(in)  :: wfn1_param
          type(wfn1_prminst), intent(in)  :: wfn1_inst
          type(wfn1_matcnst), intent(out) :: wfn1_mat
        end subroutine
        module subroutine wfn1_deallocate_mat(wfn1_mat)
          implicit none
          type(wfn1_matcnst), intent(inout) :: wfn1_mat
        end subroutine
        module subroutine wfn1_gen_mat(wfn1_inst,wfn1_mat)
          implicit none
          type(wfn1_prminst), intent(in)    :: wfn1_inst
          type(wfn1_matcnst), intent(inout) :: wfn1_mat
        end subroutine
        module subroutine wfn1_gen_ten2(wfn1_inst,wfn1_mat)
          implicit none
          type(wfn1_prminst), intent(in)    :: wfn1_inst
          type(wfn1_matcnst), intent(inout) :: wfn1_mat
        end subroutine
        module subroutine wfn1_4index_aa(wfn1_param,g_na,g_4ao,g_4mo,
     &                                   g_4t1,g_bfbf)
          implicit none
          type(wfn1_prm), intent(in) :: wfn1_param
          integer,        intent(in) :: g_na
          integer,        intent(in) :: g_4ao
          integer,        intent(in) :: g_4mo
          integer,        intent(in) :: g_4t1
          integer,        intent(in) :: g_bfbf
        end subroutine
        module subroutine wfn1_4index_ab(wfn1_param,g_na,g_nb,g_4ao,
     &                                   g_4mo,g_4t1,g_bfbf)
          implicit none
          type(wfn1_prm), intent(in) :: wfn1_param
          integer,        intent(in) :: g_na
          integer,        intent(in) :: g_nb
          integer,        intent(in) :: g_4ao
          integer,        intent(in) :: g_4mo
          integer,        intent(in) :: g_4t1
          integer,        intent(in) :: g_bfbf
        end subroutine
      end interface
      end module wfn1
c-----------------------------------------------------------------------
