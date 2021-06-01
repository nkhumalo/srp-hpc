!> \brief Scale the ranges
!>
!> Read the ranges from the command arguments, followed by scale factor,
!> scale the ranges to scale the volume by the scale factor. Output the results.
!>
      program rescale
      implicit none
      integer :: nargs, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rinp(:)
      real(kind=8),      allocatable :: rout(:)
      real(kind=8)                   :: fscale
      nargs = command_argument_count()
      if (nargs.ne.25) then
        write(*,*)"scale.x nargs = ",nargs
        stop
      endif
      allocate(argc(0:nargs))
      allocate(rinp(1:nargs))
      allocate(rout(1:nargs))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nargs
        read(argc(ii),*)rinp(ii)
      enddo
!     fscale=rinp(nargs)**(1.0d0/dble(nargs-1))
      fscale=rinp(nargs)
      do ii = 1, nargs-1
        rout(ii) = rinp(ii)*fscale
      enddo
      write(*,*)rout(1:nargs-1)
      deallocate(argc)
      deallocate(rinp)
      deallocate(rout)
      end
