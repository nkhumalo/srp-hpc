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
        module subroutine wfn1_add_occd(nmo,occ,g_occ)
          implicit none
          integer,          intent(in) :: nmo
          double precision, intent(in) :: occ(nmo)
          integer,          intent(in) :: g_occ
        end subroutine
        module subroutine wfn1_allocate_drv(wfn1_param,wfn1_deriv)
          implicit none
          type(wfn1_prm),   intent(in)  :: wfn1_param
          type(wfn1_drv),   intent(out) :: wfn1_deriv
        end subroutine
        module subroutine wfn1_allocate_mat(wfn1_param,wfn1_inst,
     &                                      wfn1_mat)
          implicit none
          type(wfn1_prm),     intent(in)  :: wfn1_param
          type(wfn1_prminst), intent(in)  :: wfn1_inst
          type(wfn1_matcnst), intent(out) :: wfn1_mat
        end subroutine
        module subroutine wfn1_allocate_wfn(wfn1_param,wfn1_inst,
     &                                      wfn1_wave)
          implicit none
          type(wfn1_prm),     intent(in)    :: wfn1_param
          type(wfn1_prminst), intent(in)    :: wfn1_inst
          type(wfn1_wfn),     intent(inout) :: wfn1_wave
        end subroutine
        module subroutine wfn1_calc_occ(g_c,nmo,ne,icnf,dcnt,occ)
          implicit none
          integer,          intent(in)  :: g_c
          integer,          intent(in)  :: nmo
          integer,          intent(in)  :: ne
          integer,          intent(in)  :: icnf(nmo)
          double precision, intent(in)  :: dcnt(nmo)
          double precision, intent(out) :: occ(nmo)
        end subroutine
        module subroutine wfn1_calc_occ2(g_c,nmo,ne,icnf,dcnt,g_occ)
          implicit none
          integer,          intent(in) :: g_c
          integer,          intent(in) :: nmo
          integer,          intent(in) :: ne
          integer,          intent(in) :: icnf(nmo)
          double precision, intent(in) :: dcnt(nmo)
          integer,          intent(in) :: g_occ
        end subroutine
        module subroutine wfn1_calc_ocx(g_c,nmo,ne,icnf,dcnt,g_ocx)
          implicit none
          integer,          intent(in) :: g_c
          integer,          intent(in) :: nmo
          integer,          intent(in) :: ne
          integer,          intent(in) :: icnf(nmo)
          double precision, intent(in) :: dcnt(nmo)
          integer,          intent(in) :: g_ocx
        end subroutine
        module subroutine wfn1_calc_docx(g_c,nmo,ne,icnf,dcnt,g_docx)
          implicit none
          integer,          intent(in) :: g_c
          integer,          intent(in) :: nmo
          integer,          intent(in) :: ne
          integer,          intent(in) :: icnf(nmo)
          double precision, intent(in) :: dcnt(nmo)
          integer,          intent(in) :: g_docx
        end subroutine
        module subroutine wfn1_configuration(wfn1_param,wfn1_wave)
          implicit none
          type(wfn1_prm), intent(in)    :: wfn1_param
          type(wfn1_wfn), intent(inout) :: wfn1_wave
        end subroutine
        module subroutine wfn1_deallocate_drv(wfn1_deriv)
          implicit none
          type(wfn1_drv),   intent(inout) :: wfn1_deriv
        end subroutine
        module subroutine wfn1_deallocate_mat(wfn1_mat)
          implicit none
          type(wfn1_matcnst), intent(inout) :: wfn1_mat
        end subroutine
        module subroutine wfn1_deallocate_wfn(wfn1_wave)
          implicit none
          type(wfn1_wfn),     intent(inout) :: wfn1_wave
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
        module subroutine wfn1_guess_corr(wfn1_param,wfn1_wave,eig,neig)
          implicit none
          type(wfn1_prm),   intent(in)    :: wfn1_param
          type(wfn1_wfn),   intent(inout) :: wfn1_wave
          integer,          intent(in)    :: neig
          double precision, intent(in)    :: eig(neig)
        end subroutine
        module subroutine wfn1_guess_corr_spin(nmo,ncore,nvirt,g_c,
     &                                         eig,neig,rndm,Tt)
          implicit none
          integer,          intent(in) :: nmo
          integer,          intent(in) :: ncore
          integer,          intent(in) :: nvirt
          integer,          intent(in) :: g_c
          integer,          intent(in) :: neig
          double precision, intent(in) :: eig(neig)
          complex(kind=8),  intent(in) :: rndm(:)
          double precision, intent(in) :: Tt
        end subroutine
        module subroutine wfn1_ortho(wfn1_mat,wfn1_wave)
          implicit none
          type(wfn1_matcnst), intent(in) :: wfn1_mat
          type(wfn1_wfn),  intent(inout) :: wfn1_wave
        end subroutine
        module subroutine wfn1_orthog(g_vecs, g_over, ometric)
          implicit none
          integer, intent(in) :: g_vecs, g_over
          logical, intent(in) :: ometric
        end subroutine
        module subroutine wfn1_orthog_dvecs(vecs, ld, ni, ulo, uhi, s,
     &                                      o, nu, w, ometric)
          implicit none
          integer, intent(in) :: ld, ni, ulo, uhi, nu
          double precision, intent(inout) ::
     &     vecs(ulo:(ulo+ld-1),1:ni), 
     &     s(ni),
     &     w(nu)
          double precision, intent(in) :: o(1:nu,ulo:(ulo+ld-1))
          logical, intent(in) :: ometric
        end subroutine
        module subroutine wfn1_orthog_zvecs(vecs, ld, ni, ulo, uhi, s,
     &                                      o, nu, w, ometric)
          implicit none
          integer, intent(in) :: ld, ni, ulo, uhi, nu
          complex(kind=8), intent(inout) ::
     &     vecs(ulo:(ulo+ld-1),1:ni), 
     &     s(ni),
     &     w(nu)
          complex(kind=8), intent(in) :: o(1:nu,ulo:(ulo+ld-1))
          logical, intent(in) :: ometric
        end subroutine
        module subroutine wfn1_tst_ortho(wfn1_mat,wfn1_wave,g_bfmo,
     &                                   gd_momo,gz_momo)
          implicit none
          type(wfn1_matcnst), intent(in) :: wfn1_mat
          type(wfn1_wfn),  intent(inout) :: wfn1_wave
          integer, intent(inout) :: g_bfmo  ! MT_DBL
          integer, intent(inout) :: gd_momo ! MT_DBL
          integer, intent(inout) :: gz_momo ! MT_DCPL
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
