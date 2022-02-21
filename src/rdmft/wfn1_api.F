!>-----------------------------------------------------------------------
!>
!> \brief Read the input
!>
!> Read the input for the RDMFT calculations
!>
      subroutine wfn1_input(rtdb)
      implicit none
      integer, intent(inout) :: rtdb !> [Input] The runtime database
      end subroutine wfn1_input
!>
!>-----------------------------------------------------------------------
!>
!> \brief Optimize the energy
!>
      function wfn1_task_energy(rtdb) result(wfn1_success)
      use wfn1
      implicit none
      integer, intent(inout) :: rtdb !> [Input] The runtime database
      logical :: wfn1_success
      wfn1_success = wfn1_task_energy_doit(rtdb)
      end function wfn1_task_energy
!>
!>-----------------------------------------------------------------------
