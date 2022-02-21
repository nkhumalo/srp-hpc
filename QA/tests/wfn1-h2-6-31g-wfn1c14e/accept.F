!> \brief Decide to accept the step (or not)
!>
!> Read the new result and the old result from the command line
!> parameters. If the new result is better than the old result
!> accept it (print "1"), else if the new result is worse accept
!> it with a certain probability. Here the probability depends
!> on how much worse the new result is compared to where we are now.
!>
      program accept
      implicit none
      integer :: nargs, ii
      character(len=32), allocatable :: argc(:)
      real(kind=8)                   :: rnew
      real(kind=8)                   :: rold
      real(kind=8)                   :: rscale
      real(kind=8) :: rrand(1)
      real(kind=8) :: p ! probability
      nargs = command_argument_count()
      if (nargs.lt.2.or.nargs.gt.3) then
        write(*,*)"accept.x nargs = ",nargs
        stop
      endif
      allocate(argc(0:nargs))
      do ii = 1, nargs
        call get_command_argument(ii,argc(ii))
      enddo
      read(argc(1),*)rnew
      read(argc(2),*)rold
      if (nargs.eq.3) then
        read(argc(3),*)rscale
      else
        rscale = 1.0d0
      endif
      if (rnew.lt.rold) then
        ! accept
        write(*,*)"1"
      else
        p = exp(-rscale*(rnew-rold)/(min(0.25d0,rold)+1.0d-6))
        call random_number(rrand)
        if (rrand(1).lt.p) then
          ! accept
          write(*,*)"1"
        else
          ! reject
          write(*,*)"0"
        endif
      endif
      deallocate(argc)
      end
