C> \ingroup wfn1
C> @{
C>
C> \brief Insert a new data point in an ordered table
C>
C> The table is sorted in increasing order. This routine adds another
C> data point to the table.
C>
      subroutine wfn1_insert(mxpts,npts,xx,xb,fx,df,xnew,bnew,
     &                       fnew,dfnew,ipos)
      implicit none
c
#include "errquit.fh"
c
!<    The dimension of the table
      integer, intent(in)    :: mxpts 
!<    The number of points in the table (must be greater than 0)
      integer, intent(inout) :: npts
!<    The coordinates \f$x_i\f$
      double precision, intent(inout) :: xx(0:mxpts)
!<    The coordinates \f$b(x_i)\f$
      double precision, intent(inout) :: xb(0:mxpts)
!<    The function values \f$f(x_i)\f$
      double precision, intent(inout) :: fx(0:mxpts)
!<    The function derivative values \f$\mathrm{d}f(x_i)\f$
      double precision, intent(inout) :: df(0:mxpts)
!<    The coordinate of the new point
      double precision, intent(in)    :: xnew
!<    The b-coordinate of the new point
      double precision, intent(in)    :: bnew
!<    The function value of the new point
      double precision, intent(in)    :: fnew
!<    The function derivative value of the new point
      double precision, intent(in)    :: dfnew
!<    The position of the new value in table
      integer, intent(out) :: ipos
c
      integer ii, jj !< Counters
c
      if (npts.ge.mxpts) call errquit("wfn1_insert: table full",0,UERR)
c
      ii = npts
      do while (xnew.lt.xx(ii).and.ii.gt.0)
        ii = ii - 1
      enddo
      if (xnew.lt.xx(ii)) then
        do jj = npts, ii, -1
          xx(jj+1) = xx(jj)
          xb(jj+1) = xb(jj)
          fx(jj+1) = fx(jj)
          df(jj+1) = df(jj)
        enddo
        xx(ii) = xnew
        xb(ii) = bnew
        fx(ii) = fnew
        df(ii) = dfnew
        ipos   = ii
      else
        do jj = npts, ii+1, -1
          xx(jj+1) = xx(jj)
          xb(jj+1) = xb(jj)
          fx(jj+1) = fx(jj)
          df(jj+1) = df(jj)
        enddo
        xx(ii+1) = xnew
        xb(ii+1) = bnew
        fx(ii+1) = fnew
        df(ii+1) = dfnew
        ipos     = ii+1
      endif
      npts = npts + 1
c
      end
C>
C> @}
