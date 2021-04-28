!> \brief Apply a random perturbation
!>
!> Read numbers from the command arguments, followed by ranges,
!> pull random numbers from the ranges and add them to the
!> input numbers. Output the results. The number of command
!> line arguments obviously has to be even.
!>
      program random
      implicit none
      integer :: nargs, nnum, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rinp(:)
      real(kind=8),      allocatable :: rrange(:)
      real(kind=8),      allocatable :: rout(:)
      real(kind=8),      allocatable :: rrand(:)
      nargs = command_argument_count()
      nnum  = nargs/2
      allocate(argc(0:nargs))
      allocate(rinp(1:nnum))
      allocate(rrange(1:nnum))
      allocate(rout(1:nnum))
      allocate(rrand(1:nnum))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nnum
        read(argc(ii),*)rinp(ii)
        read(argc(nnum+ii),*)rrange(ii)
      enddo
      call random_number(rrand)
      do ii = 1, nnum
        rrand(ii) = (2.0d0*rrand(ii)-1.0d0)*rrange(ii)
        rout(ii)  = abs(rinp(ii)+rrand(ii))
      enddo
      write(*,*)rout
      deallocate(argc)
      deallocate(rinp)
      deallocate(rrange)
      deallocate(rout)
      deallocate(rrand)
      end
