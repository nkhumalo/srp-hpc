C> \ingroup wfn1
C> @{
C>
C> \brief Fit a parabola given three points
C>
C> Given \f$ f(x_1) \f$, \f$ f(x_2) \f$, and \f$ f(x_3) \f$
C> find the parabola \f$ f(x)=ax^2+bx+c \f$ matching this
C> data, and predict its minimum. We assume that
C> \f$ x_1 \lt x_2 \lt x_3 \f$.
C>
C> The original fitting procedure has been replaced by one that is 
C> based on the lecture by Peter Young [1]. This replacement turned
C> out to be necessary because the factors \f$a\f$, \f$b\f$ and
C> \f$c\f$ can have significant errors. If these errors are not 
C> taken into account, completely false conclusions may be inferred,
C> causing the optimization to fail!
C>
C> Essentially Young's approach is based on a VanderMonde matrix.
C> It performs a least squares fit of the parabola to the available
C> data.
C>
C> ### References ###
C>
C> [1] Peter Young, "Everything you wanted to know about Data Analysis
C>     and Fitting but were afraid to ask." 2014.
C>     https://arxiv.org/pdf/1210.3781.pdf
C>
      subroutine wfn1_f1f2f3(xmin,xmax,x1,x2,x3,f1,f2,f3,a,b,c,xm,fxm)
      implicit none
c
#include "errquit.fh"
c
      double precision, intent(in)  :: xmin !< minimum \f$x\f$
      double precision, intent(in)  :: xmax !< maximum \f$x\f$
      double precision, intent(in)  :: x1  !< \f$x_1\f$
      double precision, intent(in)  :: x2  !< \f$x_2\f$
      double precision, intent(in)  :: x3  !< \f$x_3\f$
      double precision, intent(in)  :: f1  !< \f$f(x_1)\f$
      double precision, intent(in)  :: f2  !< \f$f(x_2)\f$
      double precision, intent(in)  :: f3  !< \f$f(x_3)\f$
      double precision, intent(out) :: a   !< The coefficient of \f$x^2\f$
      double precision, intent(out) :: b   !< The coefficient of \f$x\f$
      double precision, intent(out) :: c   !< The coefficient of \f$x^0\f$
      double precision, intent(out) :: xm  !< The value of \f$x\f$ that
                                           !< minimizes \f$f(x)\f$
      double precision, intent(out) :: fxm !< The value of \f$f(xm)\f$
c
      double precision xn, xu
      double precision f
      double precision aa(1:3) ! The polynomial weights 
      double precision xx(1:3) ! The positions 
      double precision ff(1:3) ! The function values
      double precision u(0:2,0:2) ! The U-matrix (Eq.(80) [1])
      double precision v(0:2)     ! The v-vector (Eq.(81) [1])
      integer, parameter :: lwork = 192
      double precision w(lwork) ! work array
      integer :: ip(3) ! pivots
      integer :: info  ! LAPACK error code
c
      integer :: ii, jj, kk ! counters
c
      double precision :: sa, sb, sc ! standard deviations
c
      f(xm) = a*xm*xm+b*xm+c
      xx(1) = x1
      xx(2) = x2
      xx(3) = x3
      ff(1) = f1
      ff(2) = f2
      ff(3) = f3
c
      u = 0.0d0
      v = 0.0d0
      do jj = 0, 2
        do ii = 0, 2
          do kk = 1, 3
            u(ii,jj) = u(ii,jj) + xx(kk)**(ii+jj)
          enddo ! kk
        enddo ! ii
        do kk = 1, 3
          v(jj) = v(jj) + ff(kk)*(xx(kk)**jj)
        enddo ! kk
      enddo ! jj
      call wfn1_dsytrf('u',3,u,3,ip,w,-1,info)
      if (w(1)-1.0d0*lwork.lt.0.0d0) then
        write(*,*)'wfn1_f1f2f3: lwork too small: ',lwork,w(1)
      endif
      call wfn1_dsytrf('u',3,u,3,ip,w,lwork,info)
      if (info.ne.0) then
        write(*,*)'wfn1_f1f2f3: dsytrf produced info = ',info
        call errquit('wfn1_f1f2f3: dsytrf failed',info,UERR)
      endif
      call wfn1_dsytri('u',3,u,3,ip,w,info)
      if (info.ne.0) then
        write(*,*)'wfn1_f1f2f3: dsytri produced info = ',info
        call errquit('wfn1_f1f2f3: dsytri failed',info,UERR)
      endif
      do jj = 0, 2
        do ii = jj+1, 2
          u(ii,jj) = u(jj,ii)
        enddo
      enddo
      call dgemv('n',3,3,1.0d0,u,3,v,1,0.0d0,aa,1)
      c  = aa(1)
      b  = aa(2)
      a  = aa(3)
      sc = dsqrt(u(0,0))
      sb = dsqrt(u(1,1))
      sa = dsqrt(u(2,2))
c
      if (a.gt.0.0d0) then
        xm = -b/(2.0d0*a)
        if (xm.lt.xmin) then
          xm = xmin
        endif
        if (xm.gt.xmax) then
          xm = xmax
        endif
      else
c
c       the function has no minimum 
c
        if (2.0d0*a*x2+b.lt.0.0d0) then
          xm = xmax
        else
          xm = xmin
        endif
      endif
      fxm = f(xm)
c
      end

C> @}
c
c-----------------------------------------------------------------------
c
#ifdef TEST_F1F2F3
      program test_wfn1_f1f2f3
      implicit none
c
c     The subroutine wfn1_f1f2f3 tries to fit a parabola through three
c     function values f(x_1), f(x_2), and f(x_3). The coordinates
c     are sorted in increasing order, i.e. x_1 <= x_2 <= x_3. In this
c     program we
c     1. Choose a parabola f(x) = a*x**2 + b*x + c by fixing a, b, and c
c     2. We choose 3 coordinates and evaluate f(x)
c     3. Call wfn1_f1f2f3
c     4. Check whether the fitted values fa, fb, and fc match the 
c        original ones.
c
      double precision :: a,  b,  c
      double precision :: da, db, dc
      double precision :: fa, fb, fc
      double precision :: x1, x2, x3
      double precision :: f1, f2, f3
      double precision :: xm, fxm
      integer :: ii
      do ii = 1, 1000
        call random_number(a)
        call random_number(b)
        call random_number(c)
        a=(a-0.5d0)*10.0d0
        b=(b-0.5d0)*10.0d0
        c=(c-0.5d0)*10.0d0
        call random_number(x1)
        call random_number(x2)
        call random_number(x3)
        x2 = x2 + x1
        x3 = x3 + x2
        x1 = x1 - 1.5d0
        x2 = x2 - 1.5d0
        x3 = x3 - 1.5d0
        f1 = a*x1**2+b*x1+c
        f2 = a*x2**2+b*x2+c
        f3 = a*x3**2+b*x3+c
        call wfn1_f1f2f3(-2.0d0,2.0d0,x1,x2,x3,f1,f2,f3,fa,fb,fc,xm,fxm)
        if (abs(a-fa).gt.1.0d-8*abs(a+fa) .or.
     &      abs(b-fb).gt.1.0d-8*abs(b+fb) .or.
     &      abs(c-fc).gt.1.0d-8*abs(c+fc)) then
          write(*,*)ii,a, b, c
          write(*,*)ii,fa,fb,fc
          write(*,*)
        endif
      enddo
      end program test_wfn1_f1f2f3
#endif
c
c-----------------------------------------------------------------------
