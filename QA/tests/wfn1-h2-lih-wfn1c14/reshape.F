!> \brief Reshape the random box
!>
!> Read the current box dimensions, the reject dimensions, and the
!> accept dimensions, at constant volume adapt the box dimensions
!> towards the accept dimensions, output the new dimensions.
!>
      program reshape_box
      implicit none
      integer :: nargs, nnum, ii, inum
      character(len=32), allocatable :: argc(:)
      real(kind=8),      allocatable :: rbox(:)
      real(kind=8),      allocatable :: rjct(:)
      real(kind=8),      allocatable :: acpt(:)
      real(kind=8),      allocatable :: rout(:)
      real(kind=8) :: volume, rescale
      nargs = command_argument_count()
      nnum  = nargs/3
      if (nargs.ne.24*3) then
        write(*,*)'reshape.x nargs = ',nargs
        stop
      endif
      allocate(argc(0:nargs))
      allocate(rbox(1:nnum))
      allocate(rjct(1:nnum))
      allocate(acpt(1:nnum))
      allocate(rout(1:nnum))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      do ii = 1, nnum
        read(argc(ii),*)rbox(ii)
        read(argc(nnum+ii),*)rjct(ii)
        read(argc(2*nnum+ii),*)acpt(ii)
      enddo
      volume=1.0d0
      inum=0
      do ii = 1, nnum
        if (rbox(ii).gt.0.0d0) then
          inum = inum + 1
          volume = volume*rbox(ii)
        endif
      enddo
      do ii = 1, nnum
        if (rjct(ii).gt.acpt(ii).and.acpt(ii).gt.0.0d0) then
          rout(ii) = rbox(ii)*sqrt(acpt(ii)/rjct(ii))
        else
          rout(ii) = rbox(ii)
        endif
      enddo
      rescale = volume
      do ii = 1, nnum
        if (rout(ii).gt.0.0d0) then
          rescale = rescale/rout(ii)
        endif
      enddo
      rescale = rescale**(1.0d0/dble(inum))
      do ii = 1, nnum
        rout(ii) = rout(ii)*rescale
      enddo
      write(*,*)rout
      deallocate(argc)
      deallocate(rbox)
      deallocate(rjct)
      deallocate(acpt)
      deallocate(rout)
      end
