!> \brief Generate a random step
!>
!> Read the ranges from the command line. Draw random 
!> numbers from these ranges. Output the results.
!>
      program random
      implicit none
      integer :: nargs, nnum, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rrange(:)
      real(kind=8),      allocatable :: rrand(:)
      nargs = command_argument_count()
      allocate(argc(0:nargs))
      allocate(rrange(1:nargs))
      allocate(rrand(1:nargs))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nnum
        read(argc(ii),*)rrange(ii)
      enddo
      call random_number(rrand)
      do ii = 1, nnum
        rrand(ii) = (2.0d0*rrand(ii)-1.0d0)*rrange(ii)
      enddo
      write(*,*)rrand
      deallocate(argc)
      deallocate(rrange)
      deallocate(rrand)
      end
