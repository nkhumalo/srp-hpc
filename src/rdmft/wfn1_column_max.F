c-----------------------------------------------------------------------
c
C> \brief Find the position of the maximum value in each column
C>
C> For a maximum overlap approach we need to know, for every column
C> vector, where it came from. Assuming that the vectors do not change
C> too much in a single iteration, the position of the absmax
C> coefficient in the column is the original position of the column
C> vector. Hence, we need to find the position of the absmax value
C> in every column.
C>
      subroutine wfn1_column_max(wfn1_param,g_a,icol)
      implicit none
!include "wfn1_param.fh"
#include "mafdecls.fh"
#include "global.fh"
#include "stdio.fh"
#include "errquit.fh"
!     Calculation parameters
      type(wfn1_prm), intent(in)     :: wfn1_param
!     Global array
      integer, intent(in) :: g_a
!     Column maximum positions
      integer, intent(out) :: icol(1:wfn1_param%nmo)
c
c     Local
c
      integer :: iproc  ! my rank
      integer :: nproc  ! total number of ranks
      integer :: nmo    ! # molecular orbitals
      integer :: ncols  ! # columns per rank
      integer :: ncollo ! the lowest column for this rank
      integer :: ncolhi ! the highest column for this rank
      integer :: ii, jj ! counters
      integer :: itype  ! type of the global array g_a
c
      logical, allocatable :: otaken(:)
      logical, allocatable :: oprior(:)
c
      double precision, allocatable :: val(:)
      complex(kind=8),  allocatable :: zval(:)
c
      character*(16) :: pname
c
      pname = "wfn1_column_max:"
      iproc = ga_nodeid()
      nproc = ga_nnodes()
      nmo   = wfn1_param%nmo
      ncols = (nmo+nproc-1)/nproc
      icol  = 0
      ncollo = iproc*ncols+1
      ncolhi = (iproc+1)*ncols
      ncollo = min(ncollo,nmo+1)
      ncolhi = min(ncolhi,nmo)
      call ga_inquire(g_a,itype,ii,jj)
      if (itype.eq.MT_DBL) then
        allocate(val(1:nmo))
      else if (itype.eq.MT_DCPL) then
        allocate(zval(1:nmo))
      else
        call errquit(pname//" invalid type",itype,UERR)
      endif
c
c     Find the position of the maximum coefficient in every column
c
      do ii = ncollo, ncolhi
        if (itype.eq.MT_DBL) then
          call ga_get(g_a,1,nmo,ii,ii,val,nmo)
        else if (itype.eq.MT_DCPL) then
          call ga_get(g_a,1,nmo,ii,ii,zval,nmo)
        endif
        icol(ii) = 1
        do jj = 2, nmo
          if (itype.eq.MT_DBL) then
            if (abs(val(jj)).gt.abs(val(icol(ii)))) then
              icol(ii) = jj
            endif
          else if (itype.eq.MT_DCPL) then
            if (abs(zval(jj)).gt.abs(zval(icol(ii)))) then
              icol(ii) = jj
            endif
          endif
        enddo
      enddo
      call ga_igop(10,icol,nmo,"+")
      if (iproc.eq.0) then
c
c       In the ideal case every maximum coefficient should come up once
c       and only once. In practice some entries might be duplicated
c       and other never used. If an entry is duplicated we need to 
c       resolve the duplication, which is done by the root process.
c
        allocate(otaken(1:nmo))
        allocate(oprior(1:nmo))
        otaken = .false.
        oprior = .false.
        do ii = 1, nmo
          otaken(icol(ii)) = .true.
        enddo
        do ii = 1, nmo
          if (oprior(icol(ii))) then
c
c           Some other column already has this as the maximum
c           coefficient. This needs to be resolved.
c
            if (itype.eq.MT_DBL) then
              call ga_get(g_a,1,nmo,ii,ii,val,nmo)
            else if (itype.eq.MT_DCPL) then
              call ga_get(g_a,1,nmo,ii,ii,zval,nmo)
            endif
            do jj = 1, nmo
              if (.not.otaken(jj)) then
                icol(ii) = jj
                exit
              endif
            enddo
            do jj = 1, nmo
              if (itype.eq.MT_DBL) then
                if (.not.otaken(jj).and.
     &              abs(val(jj)).gt.abs(val(icol(ii)))) then
                  icol(ii) = jj
                endif
              else if (itype.eq.MT_DCPL) then
                if (.not.otaken(jj).and.
     &              abs(zval(jj)).gt.abs(zval(icol(ii)))) then
                  icol(ii) = jj
                endif
              endif
            enddo
            otaken(icol(ii)) = .true.
            oprior(icol(ii)) = .true.
          else
            oprior(icol(ii)) = .true.
          endif
        enddo
        do ii = 1, nmo
          if (.not.otaken(ii)) then
            call errquit(pname//" unused coefficient ",ii)
          endif
        enddo
        deallocate(oprior)
        deallocate(otaken)
      else
        icol = 0
      endif
      if (itype.eq.MT_DBL) then
        deallocate(val)
      else if (itype.eq.MT_DCPL) then
        deallocate(zval)
      endif
      call ga_igop(10,icol,nmo,"+")
      end subroutine wfn1_column_max
c
c-----------------------------------------------------------------------
