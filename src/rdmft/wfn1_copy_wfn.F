!> \brief Copy a wavefunction
!>
      subroutine wfn1_copy_wfn(wfn_in,wfn_out)
      implicit none
      type(wfn1_wfn), intent(in)    :: wfn_in
      type(wfn1_wfn), intent(inout) :: wfn_out
      call ga_copy(wfn_in%g_na,wfn_out%g_na)
      call ga_copy(wfn_in%g_nb,wfn_out%g_nb)
      call ga_copy(wfn_in%g_ca,wfn_out%g_ca)
      call ga_copy(wfn_in%g_cb,wfn_out%g_cb)
      end subroutine wfn1_copy_wfn

