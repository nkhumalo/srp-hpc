!> \brief Absolute maximum
!>
!> Read two sets of numbers from the command arguments, 
!> take the element wise absolute maximum of the two sets,
!> output the set of maximum values.
!> The number of command line arguments obviously has to be even.
!>
      program maximum
      implicit none
      integer :: nargs, nnum, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rinp(:)
      real(kind=8),      allocatable :: rnew(:)
      real(kind=8),      allocatable :: rout(:)
      nargs = command_argument_count()
      nnum  = nargs/2
      if (nargs.ne.48) then
        write(*,*)"max.x nargs = ",nargs
        stop
      endif
      allocate(argc(0:nargs))
      allocate(rinp(1:nnum))
      allocate(rnew(1:nnum))
      allocate(rout(1:nnum))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nnum
        read(argc(ii),*)rinp(ii)
        read(argc(nnum+ii),*)rnew(ii)
      enddo
      do ii = 1, nnum
        rout(ii) = max(abs(rinp(ii)),abs(rnew(ii)))
      enddo
      write(*,*)rout
      deallocate(argc)
      deallocate(rinp)
      deallocate(rnew)
      deallocate(rout)
      end
