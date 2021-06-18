!> \brief Subtract a perturbation
!>
!> Read numbers from the command arguments, followed by a step,
!> subtract the step to the input numbers. Output the results.
!> The number of command line arguments obviously has to be even.
!>
      program subtract
      implicit none
      integer :: nargs, nnum, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rinp(:)
      real(kind=8),      allocatable :: rstep(:)
      real(kind=8),      allocatable :: rout(:)
      nargs = command_argument_count()
      nnum  = nargs/2
      allocate(argc(0:nargs))
      allocate(rinp(1:nnum))
      allocate(rstep(1:nnum))
      allocate(rout(1:nnum))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nnum
        read(argc(ii),*)rinp(ii)
        read(argc(nnum+ii),*)rstep(ii)
      enddo
      do ii = 1, nnum
        rout(ii) = rinp(ii)-rstep(ii)
      enddo
      write(*,*)rout
      deallocate(argc)
      deallocate(rinp)
      deallocate(rstep)
      deallocate(rout)
      end
