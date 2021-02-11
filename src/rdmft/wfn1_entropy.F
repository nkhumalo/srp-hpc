c-----------------------------------------------------------------------
c
      function wfn1_redist() result(fredist)
      use wfn1_fac
      implicit none
      double precision :: fredist
      fredist = fac_redist
      end function wfn1_redist
c
c-----------------------------------------------------------------------
c
C>
C> \brief Calculate the 1-electron entropy associated with an 
C> orbital and its derivative
C>
C> The entropy of a natural orbital with occupation number \f$p\f$ is
C> given by \f$S^1 = p\ln(p)+(1-p)\ln(1-p)\f$. As a result the
C> derivative is given by
C> \f$\frac{\mathrm{d}S_1}{\mathrm{d}p}=\ln(p)-\ln(1-p)\f$.
C>
C> The GSL parts are just to test a different approximation to the
C> entropy. Obviously, this introduces a licensing conflict. Therefore
C> if this expression does what I need, the solution is to expand it
C> in a Taylor series (around 0.25 as the entropy is symmetric around
C> 0.5).
C>
C> An analytic function \f$f(x)\f$ may be approximated by its Taylor
C> series as
C> \f{eqnarray}{
C> f(x) &\approx& \sum_{i=0}^{n}\frac{f^{(i)}(a)}{i!}(x-a)^i
C>      + O\left(x^{n+1}\right)
C> \f}
C> Here we are interested in the entropy given by the function
C> \f$f(x)=\ln(x!)+\ln((1-x)!)\f$ whereby \f$x \in [0,1]\f$. In order to
C> optimize the free energy we also need its gradient given by
C> \f{eqnarray}{
C> f'(x) &\approx& \sum_{i=0}^{n}\frac{f^{(i+1)}(a)}{i!}(x-a)^i
C>       + O\left(x^{n+1}\right)
C> \f}
C> Furthermore we make use of the fact that the entropy expression is
C> symmetric around \f$x=1/2\f$. This means that if we expand \f$f(x)\f$
C> around \f$a=1/4\f$ then for \f$0\le x\le 1/2\f$ we return
C> \f{eqnarray}{
C> f(x) && f'(x)
C> \f}
C> whereas for \f$1/2 \lt x \le 1\f$ we return
C> \f{eqnarray}{
C> f(1-x) && -f'(1-x)
C> \f}
C> Additionally this function is always evaluated on a finite 
C> precision computer. Hence if the Taylor series has sufficiently 
C> many terms to represent the result to machine precision then
C> the Taylor is essentially equivalent to the function itself.
C>
C> In any case the main ingredients for evaluating the Taylor series
C> of \f$f(x)\f$ and \f$f'(x)\f$ are the factors
C> \f{eqnarray}{
C> f^{(i)}(a), \; i=0,n+1
C> \f}
C> The required value of \f$n\f$ can be estimated from the form of
C> the curve suggesting that at worst
C> \f$f^{(i+1)}(x)\approx f^{(i)}(x)\f$. In addition we disregard
C> the factorial in the denominator at which point the precision 
C> achieved is determined by
C> \f{eqnarray}{
C> \max(x-a)^n < \epsilon \\\\
C> n\ln(\max(x-a)) < \ln(\epsilon) \\\\
C> n < \frac{\ln(\epsilon)}{\ln(\max(x-a))} \\\\
C> n < \frac{\ln(\epsilon)}{\ln(1/4)} \\\\
C> n \le 27
C> \f}
C>
C> In fact, instead of introducing a Taylor series there is a much
C> simpler way to approximate the entropy expression
C> \f$f(x)=\ln(x!)+\ln((1-x)!)\f$ as 
C> \f$f(x)\approx [x(1-x)]^b,\; b=0.9757\f$
C> 
      subroutine wfn1_entropy1(p,s,ds)
      implicit none
!>    The occupation number
      double precision, intent(in)  :: p
!>    The entropy
      double precision, intent(out) :: s
!>    The entropy gradient
      double precision, intent(out) :: ds
c
      double precision :: a0, a
      double precision :: tol
      double precision :: ss, dss
      integer, parameter :: mxcc = 28
      double precision, parameter :: cc(0:mxcc) = (/
     &  -1.8267295744229717e-01, !   0
     &  -4.7492598692312660e-01, !   1
     &   1.9614310244009392e+00, !   2
     &  -7.6584746194002107e-01, !   3
     &   3.5824935810760326e+00, !   4
     &  -6.7105852916520634e+00, !   5
     &   3.7028427428299736e+01, !   6
     &  -1.3867070561680774e+02, !   7
     &   9.1270050624592193e+02, !   8
     &  -5.1733200953013929e+03, !   9
     &   4.0438498887562106e+04, !  10
     &  -3.0445288613979606e+05, !  11
     &   2.7940653844185043e+06, !  12
     &  -2.6013424905976593e+07, !  13
     &   2.7640984169809264e+08, !  14
     &  -3.0480304585382962e+09, !  15
     &   3.6979925453039085e+10, !  16
     &  -4.6961506675232050e+11, !  17
     &   6.4226726567984150e+12, !  18
     &  -9.2114814596364359e+13, !  19
     &   1.4041599596566502e+15, !  20
     &  -2.2420499819440248e+16, !  21
     &   3.7721542091731910e+17, !  22
     &  -6.6320480845454735e+18, !  23
     &   1.2212077163997869e+20, !  24
     &  -2.3434678559572684e+21, !  25
     &   4.6887208088579196e+22, !  26
     &  -9.7498848636924539e+23, !  27
     &   2.1063843942817377e+25/)!  28
      integer :: ii ! counter
      double precision :: x
c
#ifdef USE_GSL
      double precision, external :: gsl_sf_lngamma
      double precision, external :: gsl_sf_psi
c     a  = log(0.5d0)/(2.0d0*gsl_sf_lngamma(1.5d0))
      a  = 0.2869408590744988D+01
      s  = a*(gsl_sf_lngamma(1.0d0+p) + gsl_sf_lngamma(2.0d0-p))
      ds = a*(gsl_sf_psi(1.0d0+p) - gsl_sf_psi(2.0d0-p))
c
      ss  = 0.0d0
      dss = 0.0d0
      a0  = 0.25d0
      if (p.le.0.5d0) then
        x = p
      else
        x = 1.0d0-p
      endif
      do ii = mxcc-1, 0, -1
        ss  = ss/dble(ii+1)*(x-a0) +cc(ii)
        dss = dss/dble(ii+1)*(x-a0)+cc(ii+1)
      enddo
      if (p.le.0.5d0) then
      else
        dss = -dss
      endif
      ss  = a*ss
      dss = a*dss
c     write(*,*)"HVD s ",s, abs(s-ss),
c    &          2.0d0*abs(s-ss)/(abs(s)+abs(ss))
c     write(*,*)"HVD ds",ds,abs(ds-dss),
c    &          2.0d0*abs(ds-dss)/(abs(ds)+abs(dss))
c     write(*,'(" HVD a ",e24.16)')a
      s  = ss
      ds = dss
#elif USE_MERMIN
      tol = 1.0d-10
      if   (p.lt.tol) then
        s  = 0.0d0
        ds = log(tol)
      else if (p.gt.1.0d0-tol) then
        s  = 0.0d0
        ds = -log(tol)
      else
        s  = p*log(p) + (1.0d0-p)*log(1.0d0-p)
        ds = log(p) - log(1.0d0-p)
      endif
#else
      double precision :: b = 0.9757d0
      a = (0.5d0*log(0.5d0)+(1.0d0-0.5d0)*log(1.0d0-0.5d0))/
     &    ((0.5d0*(1.0d0-0.5d0))**b)
      tol = 1.0d-4
      if   (p.lt.tol) then
        s  = 0.0d0
        ds = a*b*((tol*(1.0d0-tol))**(b-1.0d0))*(1.0d0-2.0d0*tol)
      else if (p.gt.1.0d0-tol) then
        s  = 0.0d0
        ds = -a*b*((tol*(1.0d0-tol))**(b-1.0d0))*(1.0d0-2.0d0*tol)
      else
        s  = a*((p*(1.0d0-p))**b)
        ds = a*b*((p*(1.0d0-p))**(b-1.0d0))*(1.0d0-2.0d0*p)
      endif
#endif
c
      end subroutine wfn1_entropy1
c
c-----------------------------------------------------------------------
c
C>
C> \brief Evaluate the 1-electron entropy and its gradient but without
C> the Stirling approximation.
C>
C> The entropy of a natural orbital with occupation number \f$p\f$ is
C> given by
C> \f{eqnarray}{
C>   S^1(p) &=& \ln(p!)+\ln([1-p]!).
C> \f}
C> This equation can be expressed in terms of the Gamma-function as
C> \f{eqnarray}{
C>   S^1(p) &=& \ln(\Gamma(1+p))+\ln(\Gamma(1+[1-p]))
C> \f}
C> where \f$\ln(\Gamma(x))\f$ is also know as the log-gamma function.
C> The first order derivative of the log-gamma function is also known
C> as the digamma function, denoted by \f$\psi_0(x)\f$. In terms of the
C> digamma function the first order derivative of \f$S^1(p)\f$ is
C> \f{eqnarray}{
C>   \frac{\partial S_1}{\partial p}
C>   &=& \psi_0(1+p)-\psi_0(1+[1-p])
C> \f}
C>
C> This form of the 1-electron entropy is an essential ingredient of
C> the electron-pair entropy for the correlation model. This as opposed
C> to the finite temperature 1-electron entropy.
C>
      subroutine wfn1_entropy_si(pow,p,s,ds)
      implicit none
!>    The power in the entropy expression (akin to Renyi)
      double precision, intent(in)  :: pow
!>    The occupation number
      double precision, intent(in)  :: p
!>    The entropy
      double precision, intent(out) :: s
!>    The entropy gradient
      double precision, intent(out) :: ds
c
      double precision :: a0, a
      double precision :: tol
      double precision :: ss, dss
      integer, parameter :: mxcc = 28
      double precision, parameter :: cc(0:mxcc) = (/
     &  -1.8267295744229717e-01, !   0
     &  -4.7492598692312660e-01, !   1
     &   1.9614310244009392e+00, !   2
     &  -7.6584746194002107e-01, !   3
     &   3.5824935810760326e+00, !   4
     &  -6.7105852916520634e+00, !   5
     &   3.7028427428299736e+01, !   6
     &  -1.3867070561680774e+02, !   7
     &   9.1270050624592193e+02, !   8
     &  -5.1733200953013929e+03, !   9
     &   4.0438498887562106e+04, !  10
     &  -3.0445288613979606e+05, !  11
     &   2.7940653844185043e+06, !  12
     &  -2.6013424905976593e+07, !  13
     &   2.7640984169809264e+08, !  14
     &  -3.0480304585382962e+09, !  15
     &   3.6979925453039085e+10, !  16
     &  -4.6961506675232050e+11, !  17
     &   6.4226726567984150e+12, !  18
     &  -9.2114814596364359e+13, !  19
     &   1.4041599596566502e+15, !  20
     &  -2.2420499819440248e+16, !  21
     &   3.7721542091731910e+17, !  22
     &  -6.6320480845454735e+18, !  23
     &   1.2212077163997869e+20, !  24
     &  -2.3434678559572684e+21, !  25
     &   4.6887208088579196e+22, !  26
     &  -9.7498848636924539e+23, !  27
     &   2.1063843942817377e+25/)!  28
      integer :: ii ! counter
      double precision :: x
c
      double precision, external :: gsl_sf_lngamma
      double precision, external :: gsl_sf_psi
c     a  = log(0.5d0)/(2.0d0*gsl_sf_lngamma(1.5d0))
c     a  = 0.2869408590744988D+01
      a  = 1.0d0
#ifdef LOGGAMMA
      s  = a*(gsl_sf_lngamma(1.0d0+p) + gsl_sf_lngamma(2.0d0-p))
      ds = a*(gsl_sf_psi(1.0d0+p) - gsl_sf_psi(2.0d0-p))
#else
      s  = 0.0d0
      ds = 0.0d0
      if (pow.ge.1.0d0.and.p.gt.1.0d-12) then
        s  = (p**pow)*log(p**pow)
        ds = pow*(p**(pow-1.0d0))*log(p**pow)+pow*(p**(pow-1.0d0))
      endif
#endif
c     s  = a*(gsl_sf_lngamma(1.0d0+p))
c     ds = a*(gsl_sf_psi(1.0d0+p))
c
c     ss  = 0.0d0
c     dss = 0.0d0
c     a0  = 0.25d0
c     if (p.le.0.5d0) then
c       x = p
c     else
c       x = 1.0d0-p
c     endif
c     do ii = mxcc-1, 0, -1
c       ss  = ss/dble(ii+1)*(x-a0) +cc(ii)
c       dss = dss/dble(ii+1)*(x-a0)+cc(ii+1)
c     enddo
c     if (p.le.0.5d0) then
c     else
c       dss = -dss
c     endif
c     ss  = a*ss
c     dss = a*dss
c     write(*,*)"HVD s ",s, abs(s-ss),
c    &          2.0d0*abs(s-ss)/(abs(s)+abs(ss))
c     write(*,*)"HVD ds",ds,abs(ds-dss),
c    &          2.0d0*abs(ds-dss)/(abs(ds)+abs(dss))
c     write(*,'(" HVD a ",e24.16)')a
c     s  = ss
c     ds = dss
c
      end subroutine wfn1_entropy_si
c
c-----------------------------------------------------------------------
c
C>
C> \brief Evaluate the 2-electron entropy and its gradient but without
C> the Stirling approximation.
C>
C> This form of the 2-electron entropy is an essential ingredient of
C> the electron-pair entropy for the correlation model. 
C>
C> The entropy of a natural orbital pair with occupation numbers
C> \f$p_i\f$ and \f$p_j\f$ is
C> given by
C> \f{eqnarray}{
C>    \f$S^2 &=& \ln(p_i p_j!) 
C>               +\ln([1-p_i]p_j!) \\\\
C>           &&  +\ln(p_i[1-p_j]!)
C>               +\ln([1-p_i][1-p_j]!)
C> \f}.
C> This expression can be rewritten in terms of the log-gamma function
C> as
C> \f{eqnarray}{
C>    \f$S^2 &=& \ln(\Gamma(1+p_i p_j)) 
C>               +\ln(\Gamma(1+[1-p_i]p_j)) \\\\
C>           &&  +\ln(\Gamma(1+p_i[1-p_j]))
C>               +\ln(\Gamma(1+[1-p_i][1-p_j]))
C> \f}.
C> The derivative of \f$\ln(\Gamma(x))\f$ is know as the digamma
C> function denoted as \f$\psi_0(x)\f$.
C>
      subroutine wfn1_entropy_sij(pow,pi,pj,s,dsi,dsj)
      implicit none
!>    The power in the entropy (akin to Renyi)
      double precision, intent(in)  :: pow
!>    The occupation number for natural orbital \f$i\f$
      double precision, intent(in)  :: pi
!>    The occupation number for natural orbital \f$j\f$
      double precision, intent(in)  :: pj
!>    The entropy
      double precision, intent(out) :: s
!>    The entropy gradient wrt \f$p_i\f$
      double precision, intent(out) :: dsi
!>    The entropy gradient wrt \f$p_j\f$
      double precision, intent(out) :: dsj
c
      double precision :: a0, a
      double precision :: tol
      double precision :: ss, dss
      integer, parameter :: mxcc = 28
      double precision, parameter :: cc(0:mxcc) = (/
     &  -1.8267295744229717e-01, !   0
     &  -4.7492598692312660e-01, !   1
     &   1.9614310244009392e+00, !   2
     &  -7.6584746194002107e-01, !   3
     &   3.5824935810760326e+00, !   4
     &  -6.7105852916520634e+00, !   5
     &   3.7028427428299736e+01, !   6
     &  -1.3867070561680774e+02, !   7
     &   9.1270050624592193e+02, !   8
     &  -5.1733200953013929e+03, !   9
     &   4.0438498887562106e+04, !  10
     &  -3.0445288613979606e+05, !  11
     &   2.7940653844185043e+06, !  12
     &  -2.6013424905976593e+07, !  13
     &   2.7640984169809264e+08, !  14
     &  -3.0480304585382962e+09, !  15
     &   3.6979925453039085e+10, !  16
     &  -4.6961506675232050e+11, !  17
     &   6.4226726567984150e+12, !  18
     &  -9.2114814596364359e+13, !  19
     &   1.4041599596566502e+15, !  20
     &  -2.2420499819440248e+16, !  21
     &   3.7721542091731910e+17, !  22
     &  -6.6320480845454735e+18, !  23
     &   1.2212077163997869e+20, !  24
     &  -2.3434678559572684e+21, !  25
     &   4.6887208088579196e+22, !  26
     &  -9.7498848636924539e+23, !  27
     &   2.1063843942817377e+25/)!  28
      integer :: ii ! counter
      double precision :: x
c
      double precision, external :: gsl_sf_lngamma
      double precision, external :: gsl_sf_psi
#ifdef LOGGAMMA
      s   = gsl_sf_lngamma(1.0d0+pi*pj)
     &    + gsl_sf_lngamma(1.0d0+(1.0d0-pi)*pj)
     &    + gsl_sf_lngamma(1.0d0+(1.0d0-pj)*pi)
     &    + gsl_sf_lngamma(1.0d0+(1.0d0-pi)*(1.0d0-pj))
      dsi = gsl_sf_psi(1.0d0+pi*pj)*pj
     &    - gsl_sf_psi(1.0d0+(1.0d0-pi)*pj)*pj
     &    + gsl_sf_psi(1.0d0+(1.0d0-pj)*pi)*(1.0d0-pj)
     &    - gsl_sf_psi(1.0d0+(1.0d0-pi)*(1.0d0-pj))*(1.0d0-pj)
      dsj = gsl_sf_psi(1.0d0+pi*pj)*pi
     &    + gsl_sf_psi(1.0d0+(1.0d0-pi)*pj)*(1.0d0-pi)
     &    - gsl_sf_psi(1.0d0+(1.0d0-pj)*pi)*pi
     &    - gsl_sf_psi(1.0d0+(1.0d0-pi)*(1.0d0-pj))*(1.0d0-pi)
#else
c     von Neumann entropy
      s   = 0.0d0
      dsi = 0.0d0
      dsj = 0.0d0
      if (pow.ge.1.0d0.or.pi*pj.gt.1.0d-12) then
        s   = ((pi*pj)**pow)*log((pi*pj)**pow)
        dsi = pow*pj*((pi*pj)**(pow-1.0d0))*log((pi*pj)**pow)
     &      + pow*pj*((pi*pj)**(pow-1.0d0))
        dsj = pow*pi*((pi*pj)**(pow-1.0d0))*log((pi*pj)**pow)
     &      + pow*pi*((pi*pj)**(pow-1.0d0))
      endif
#endif
c
c     ss  = 0.0d0
c     dss = 0.0d0
c     a0  = 0.25d0
c     if (p.le.0.5d0) then
c       x = p
c     else
c       x = 1.0d0-p
c     endif
c     do ii = mxcc-1, 0, -1
c       ss  = ss/dble(ii+1)*(x-a0) +cc(ii)
c       dss = dss/dble(ii+1)*(x-a0)+cc(ii+1)
c     enddo
c     if (p.le.0.5d0) then
c     else
c       dss = -dss
c     endif
c     ss  = a*ss
c     dss = a*dss
c     write(*,*)"HVD s ",s, abs(s-ss),
c    &          2.0d0*abs(s-ss)/(abs(s)+abs(ss))
c     write(*,*)"HVD ds",ds,abs(ds-dss),
c    &          2.0d0*abs(ds-dss)/(abs(ds)+abs(dss))
c     write(*,'(" HVD a ",e24.16)')a
c     s  = ss
c     ds = dss
c
      end subroutine wfn1_entropy_sij
c
c-----------------------------------------------------------------------
c
!> \brief Approximation to the entanglement entropy
!>
!> Rather than implementing the entanglement by combining the routines
!> wfn1_entropy_si and wfn1_entropy_sij we use an approximation 
!> similar to wfn1_entropy1. The approximation is generalized to all
!> of the entanglement entropy.
!> 
      subroutine wfn1_entropy2_sij(pi,pj,s,dsi,dsj)
      implicit none
!>    The occupation number for natural orbital \f$i\f$
      double precision, intent(in)  :: pi
!>    The occupation number for natural orbital \f$j\f$
      double precision, intent(in)  :: pj
!>    The entropy
      double precision, intent(out) :: s
!>    The entropy gradient wrt \f$p_i\f$
      double precision, intent(out) :: dsi
!>    The entropy gradient wrt \f$p_j\f$
      double precision, intent(out) :: dsj
c
      double precision :: a
      double precision :: tol
      double precision :: di, dj ! local versions of pi and pj
c     double precision, parameter :: b = 0.9513d0
      double precision, parameter :: b = 0.0500d0
c
      a =  1.0d0/((0.5d0*(1.0d0-0.5d0))**(2.0d0*b))
c
      di  = max(0.0d0,pi)
      dj  = max(0.0d0,pj)
      s   = (di*(1.0d0-di)*dj*(1.0d0-dj))**b
      di  = max(1.0d-10,pi)
      dj  = max(1.0d-10,pj)
      dsi = b*(di*(1.0d0-di)*dj*(1.0d0-dj))**(b-1.0d0)
     &      *(1.0d0-2.0d0*di)*dj*(1.0d0-dj)
      dsj = b*(di*(1.0d0-di)*dj*(1.0d0-dj))**(b-1.0d0)
     &      *di*(1.0d0-di)*(1.0d0-2.0d0*dj)
c
      s   = a*s
      dsi = a*dsi
      dsj = a*dsj
c
      end subroutine wfn1_entropy2_sij
c
c-----------------------------------------------------------------------
c
!> \brief Set the parameters that define the entropy
!>
!> This subroutine is a temporary feature until there is more clarity
!> on the required parameters.
!>
      subroutine wfn1_set_entropy_fac(faa,fab,taa,dab,
     &                                px1aa,px2aa,px1ab,px2ab,
     &                                py1aa,py1ab,py2aa,py2ab,
     &                                py3aa,py3ab,py4aa,py4ab,
     &                                py5aa,py5ab,
     &                                pz1aa,pz1ab,pz2aa,pz2ab,
     &                                pz3aa,pz3ab,pz4aa,pz4ab,
     &                                pz5aa,pz5ab,
     &                                hy1aa,hy1ab,hy2aa,hy2ab,
     &                                hy3aa,hy3ab,hy4aa,hy4ab,
     &                                hy5aa,hy5ab,
     &                                hz1aa,hz1ab,hz2aa,hz2ab,
     &                                hz3aa,hz3ab,hz4aa,hz4ab,
     &                                hz5aa,hz5ab)
      use wfn1_fac
      implicit none
!>    The same-spin overal scale factor
      double precision, intent(in) :: faa
!>    The opposite-spin overal scale factor
      double precision, intent(in) :: fab
!>    The same-spin coordinate offset term
      double precision, intent(in) :: taa
!>    The opposite-spin same orbital scale factor
      double precision, intent(in) :: dab
!>    The same-spin power-x
      double precision, intent(in) :: px1aa
      double precision, intent(in) :: px2aa
!>    The opposite-spin power-x
      double precision, intent(in) :: px1ab
      double precision, intent(in) :: px2ab
!>    The same-spin power-y
      double precision, intent(in) :: py1aa
      double precision, intent(in) :: py2aa
      double precision, intent(in) :: py3aa
      double precision, intent(in) :: py4aa
      double precision, intent(in) :: py5aa
!>    The opposite-spin power-y
      double precision, intent(in) :: py1ab
      double precision, intent(in) :: py2ab
      double precision, intent(in) :: py3ab
      double precision, intent(in) :: py4ab
      double precision, intent(in) :: py5ab
!>    The same-spin power-z
      double precision, intent(in) :: pz1aa
      double precision, intent(in) :: pz2aa
      double precision, intent(in) :: pz3aa
      double precision, intent(in) :: pz4aa
      double precision, intent(in) :: pz5aa
!>    The opposite-spin power-z
      double precision, intent(in) :: pz1ab
      double precision, intent(in) :: pz2ab
      double precision, intent(in) :: pz3ab
      double precision, intent(in) :: pz4ab
      double precision, intent(in) :: pz5ab
!>    The same-spin hf power-y
      double precision, intent(in) :: hy1aa
      double precision, intent(in) :: hy2aa
      double precision, intent(in) :: hy3aa
      double precision, intent(in) :: hy4aa
      double precision, intent(in) :: hy5aa
!>    The opposite-spin hf power-y
      double precision, intent(in) :: hy1ab
      double precision, intent(in) :: hy2ab
      double precision, intent(in) :: hy3ab
      double precision, intent(in) :: hy4ab
      double precision, intent(in) :: hy5ab
!>    The same-spin hf power-z
      double precision, intent(in) :: hz1aa
      double precision, intent(in) :: hz2aa
      double precision, intent(in) :: hz3aa
      double precision, intent(in) :: hz4aa
      double precision, intent(in) :: hz5aa
!>    The opposite-spin hf power-z
      double precision, intent(in) :: hz1ab
      double precision, intent(in) :: hz2ab
      double precision, intent(in) :: hz3ab
      double precision, intent(in) :: hz4ab
      double precision, intent(in) :: hz5ab
!>
c     double precision :: fac_aa
c     double precision :: fac_ab
c     double precision :: trm_aa
c     double precision :: trm_ab
c     double precision :: pow_y_aa
c     double precision :: pow_y_ab
c     double precision :: pow_z_aa
c     double precision :: pow_z_ab
c     double precision, allocatable :: orb_ea(:)
c     double precision, allocatable :: orb_eb(:)
c     common/wfn1_fac/fac_aa,fac_ab,
c    +                trm_aa,trm_ab,
c    +                pow_y_aa,pow_y_ab,
c    +                pow_z_aa,pow_z_ab,
c    +                orb_ea,orb_eb
      fac_aa    = faa
      fac_ab    = fab
      fac_redist= taa
      del_ab    = dab
      pow_x1_aa = px1aa
      pow_x2_aa = px2aa
      pow_x1_ab = px1ab
      pow_x2_ab = px2ab
c
      pow_y1_aa = py1aa
      pow_y1_ab = py1ab
      pow_y2_aa = py2aa
      pow_y2_ab = py2ab
      pow_y3_aa = py3aa
      pow_y3_ab = py3ab
      pow_y4_aa = py4aa
      pow_y4_ab = py4ab
      pow_y5_aa = py5aa
      pow_y5_ab = py5ab
      pow_z1_aa = pz1aa
      pow_z1_ab = pz1ab
      pow_z2_aa = pz2aa
      pow_z2_ab = pz2ab
      pow_z3_aa = pz3aa
      pow_z3_ab = pz3ab
      pow_z4_aa = pz4aa
      pow_z4_ab = pz4ab
      pow_z5_aa = pz5aa
      pow_z5_ab = pz5ab
c
      powhf_y1_aa = hy1aa
      powhf_y1_ab = hy1ab
      powhf_y2_aa = hy2aa
      powhf_y2_ab = hy2ab
      powhf_y3_aa = hy3aa
      powhf_y3_ab = hy3ab
      powhf_y4_aa = hy4aa
      powhf_y4_ab = hy4ab
      powhf_y5_aa = hy5aa
      powhf_y5_ab = hy5ab
      powhf_z1_aa = hz1aa
      powhf_z1_ab = hz1ab
      powhf_z2_aa = hz2aa
      powhf_z2_ab = hz2ab
      powhf_z3_aa = hz3aa
      powhf_z3_ab = hz3ab
      powhf_z4_aa = hz4aa
      powhf_z4_ab = hz4ab
      powhf_z5_aa = hz5aa
      powhf_z5_ab = hz5ab
      end subroutine wfn1_set_entropy_fac
c
c-----------------------------------------------------------------------
c
!> \brief Set the orbital energies that define the entropy
!>
!> This subroutine is a temporary feature until there is more clarity
!> on the required parameters.
!>
      subroutine wfn1_set_entropy_eorb(eorb_a,eorb_b,nmo)
      use wfn1_fac
      implicit none
!>    The number of MOs
      integer, intent(in)          :: nmo
!>    The alpha-spin orbital energies
      double precision, intent(in) :: eorb_a(1:nmo)
!>    The beta-spin orbital energies
      double precision, intent(in) :: eorb_b(1:nmo)
!>
c     double precision :: fac_aa
c     double precision :: fac_ab
c     double precision :: trm_aa
c     double precision :: trm_ab
c     double precision :: pow_y_aa
c     double precision :: pow_y_ab
c     double precision :: pow_z_aa
c     double precision :: pow_z_ab
c     double precision, allocatable :: orb_ea(:)
c     double precision, allocatable :: orb_eb(:)
c     common/wfn1_fac/fac_aa,fac_ab,
c    +                trm_aa,trm_ab,
c    +                pow_y_aa,pow_y_ab,
c    +                pow_z_aa,pow_z_ab,
c    +                orb_ea,orb_eb
      if (.not.allocated(orb_ea)) then
        allocate(orb_ea(1:nmo))
      endif
      if (.not.allocated(orb_eb)) then
        allocate(orb_eb(1:nmo))
      endif
      orb_ea = eorb_a
      orb_eb = eorb_b
cDEBUG
      write(110,*)abs(eorb_a(1)-eorb_a(2))
      write(111,*)abs(eorb_a(1)-eorb_a(2))
      close(110)
      close(111)
cDEBUG
      end subroutine wfn1_set_entropy_eorb
c
c-----------------------------------------------------------------------
c
!> \brief Set the Hartree-Fock orbital energies that define the entropy
!>
!> This subroutine is a temporary feature until there is more clarity
!> on the required parameters.
!>
      subroutine wfn1_set_entropy_eorb_hf(eorb_a,eorb_b,nmo)
      use wfn1_fac
      implicit none
!>    The number of MOs
      integer, intent(in)          :: nmo
!>    The alpha-spin orbital energies
      double precision, intent(in) :: eorb_a(1:nmo)
!>    The beta-spin orbital energies
      double precision, intent(in) :: eorb_b(1:nmo)
!>
c     double precision :: fac_aa
c     double precision :: fac_ab
c     double precision :: trm_aa
c     double precision :: trm_ab
c     double precision :: pow_y_aa
c     double precision :: pow_y_ab
c     double precision :: pow_z_aa
c     double precision :: pow_z_ab
c     double precision, allocatable :: orb_ea(:)
c     double precision, allocatable :: orb_eb(:)
c     common/wfn1_fac/fac_aa,fac_ab,
c    +                trm_aa,trm_ab,
c    +                pow_y_aa,pow_y_ab,
c    +                pow_z_aa,pow_z_ab,
c    +                orb_ea,orb_eb
      if (.not.allocated(orb_ea_hf)) then
        allocate(orb_ea_hf(1:nmo))
      endif
      if (.not.allocated(orb_eb_hf)) then
        allocate(orb_eb_hf(1:nmo))
      endif
      orb_ea_hf = eorb_a
      orb_eb_hf = eorb_b
cDEBUG
      write(120,*)abs(eorb_a(1)-eorb_a(2))
      write(121,*)abs(eorb_a(1)-eorb_a(2))
      close(120)
      close(121)
cDEBUG
      end subroutine wfn1_set_entropy_eorb_hf
c
c-----------------------------------------------------------------------
c
C>
C> \brief Calculate the 2-electron entropy associated with a 
C> natural orbital pair and its derivatives
C>
C> The entropy for a pair of natural orbitals is given by
C> \f{eqnarray}{
C>   S(p_i,p_j) &=& S_1(p_i) + S_1(p_j) - S_2(p_i,p_j)
C> \f}
C> The normalization of this expression should be such that 
C> \f$S(1/2,1/2) = 1\f$.
C> 
      subroutine wfn1_entropy2_ab(ii,jj,pi,pj,s,dsi,dsj)
      use wfn1_fac
      implicit none
C>    The orbital index \f$i\f$
      integer, intent(in) :: ii
C>    The orbital index \f$j\f$
      integer, intent(in) :: jj
C>    The occupation number of orbital \f$i\f$
      double precision, intent(in)  :: pi
C>    The occupation number of orbital \f$i\f$
      double precision, intent(in)  :: pj
C>    The entropy
      double precision, intent(out) :: s
C>    The derivative wrt \f$p_i\f$
      double precision, intent(out) :: dsi
C>    The derivative wrt \f$p_j\f$
      double precision, intent(out) :: dsj
c
c     Local
c
      double precision :: a, ai, aj, aij
      double precision :: d, si, sj, sij
      double precision :: pii, pjj
      double precision :: ds1_i
      double precision :: ds1_j
      double precision :: ds2_i
      double precision :: ds2_j
      double precision :: f   ! the enhancement factor
      double precision :: dfi ! the gradient of the enhancement factor
      double precision :: dfj ! the gradient of the enhancement factor

c
c     Global
c
c     double precision :: fac_aa
c     double precision :: fac_ab
c     double precision :: trm_aa
c     double precision :: trm_ab
c     double precision :: pow_y_aa
c     double precision :: pow_y_ab
c     double precision :: pow_z_aa
c     double precision :: pow_z_ab
c     double precision, allocatable :: orb_ea(:)
c     double precision, allocatable :: orb_eb(:)
c     common/wfn1_fac/fac_aa,fac_ab,
c    +                trm_aa,trm_ab,
c    +                pow_y_aa,pow_y_ab,
c    +                pow_z_aa,pow_z_ab,
c    +                orb_ea,orb_eb
      double precision :: de
      double precision :: de_hf
c
#ifdef OLD
c     call wfn1_entropy_si(0.5d0,ai,ds1_i)
c     call wfn1_entropy_si(0.5d0,aj,ds1_j)
c     call wfn1_entropy_sij(0.5d0,0.5d0,aij,ds2_i,ds2_j)
#ifdef LOGGAMMA
      call wfn1_entropy_si(1.0d0/3.0d0,ai,ds1_i)
      call wfn1_entropy_si(1.0d0/3.0d0,aj,ds1_j)
      call wfn1_entropy_sij(1.0d0/3.0d0,1.0d0/3.0d0,aij,ds2_i,ds2_j)
#else
      call wfn1_entropy_si(pow_ab,0.5d0,ai,ds1_i)
      call wfn1_entropy_si(pow_ab,0.5d0,aj,ds1_j)
      call wfn1_entropy_sij(pow_ab,0.5d0,0.5d0,aij,ds2_i,ds2_j)
#endif
      a = 0.25d0/(ai+aj-aij)
c
c     The following factor comes from a H2 calculation in the STO-3G
c     basis set at an interatomic distance of the equilibrium bond-
c     length plus 5 Angstrom, and then dividing the UHF 2-electron
c     energy by the RHF 2-electron energy and subtracting that fraction
c     (0.21364) from 1. I.e. electron correlation does not annihilate
c     the 2-electron interaction completely but mostly.
c
c     a = a*0.78636d0
c     if (pi.lt.0.5d0.and.pj.lt.0.5d0) then
c       a = fac_ab_ll*a
c     else if (pi.lt.0.5d0.and.pj.ge.0.5d0) then
c       a = fac_ab_hl*a
c     else if (pi.ge.0.5d0.and.pj.lt.0.5d0) then
c       a = fac_ab_hl*a
c     else
c       a = fac_ab_hh*a
c     endif
      a = fac_ab*a
c
      call wfn1_entropy_si(pow_ab,pi,si,ds1_i)
      call wfn1_entropy_si(pow_ab,pj,sj,ds1_j)
      call wfn1_entropy_sij(pow_ab,pi,pj,sij,ds2_i,ds2_j)
c
      s   = a*(si+sj-sij)
cDEBUG
c     if (s.gt.1.0d0) then
c       write(*,*)'S too large: ',s,pi,pj
c     else if (s.lt.0.0d0) then
c       write(*,*)'S too small: ',s,pi,pj
c     endif
cDEBUG
      dsi = a*(ds1_i-ds2_i)
      dsj = a*(ds1_j-ds2_j)
c
#else
c     call wfn1_entropy2_sij(0.5d0,0.5d0,sij,dsi,dsj)
c     write(*,*)'HVD: sij = ',sij,fac_ab
      a = fac_ab
      d = del_ab
cDEBUG
c     write(*,*)'HVD: fac_ab = ',fac_ab,a
c     write(*,*)'HVD: pi,pj  = ',pi,pj
cDEBUG
      if      (pi.gt.0.5d0.and.pj.gt.0.5d0) then
        f   = 1.0d0
        dfi = 0.0d0
        dfj = 0.0d0
      else if (pi.lt.0.5d0.and.pj.lt.0.5d0) then
c       call wfn1_entropy2_ll(pi,pj,f,dfi,dfj)
        f   = 1.0d0
        dfi = 0.0d0
        dfj = 0.0d0
      else
        f   = 1.0d0
        dfi = 0.0d0
        dfj = 0.0d0
      endif
c     call wfn1_entropy2_sij(pi,pj,s,dsi,dsj)
      if (pi.lt.1.0d-10) then
        pii = 1.0d-10
      else if (pi.gt.1.0d0-1.0d-10) then
        pii = 1.0d0-1.0d-10
      else
        pii = pi
      endif
      if (pj.lt.1.0d-10) then
        pjj = 1.0d-10
      else if (pj.gt.1.0d0-1.0d-10) then
        pjj = 1.0d0-1.0d-10
      else
        pjj = pj
      endif
cDEBUG
c     write(0,*)'HVD: ii,jj,Ea_ii,Eb_jj=',ii,jj,orb_ea(ii),orb_eb(jj)
cDEBUG
      de    = orb_ea(ii)-orb_eb(jj)
      de_hf = orb_ea_hf(ii)-orb_eb_hf(jj)
      call wfn1_thimble(pii,pjj,de_hf,de,
     &                  pow_x1_ab,pow_x2_ab,
     &                  pow_y1_ab,pow_y2_ab,
     &                  pow_y3_ab,pow_y4_ab,
     &                  pow_y5_ab,
     &                  pow_z1_ab,pow_z2_ab,
     &                  pow_z3_ab,pow_z4_ab,
     &                  pow_z5_ab,
     &                  powhf_y1_ab,powhf_y2_ab,
     &                  powhf_y3_ab,powhf_y4_ab,
     &                  powhf_y5_ab,
     &                  powhf_z1_ab,powhf_z2_ab,
     &                  powhf_z3_ab,powhf_z4_ab,
     &                  powhf_z5_ab,
     &                  s,dsi,dsj)
      if (ii.eq.jj) then
        s   = d*s*f
        dsi = d*(f*dsi+s*dfi)
        dsj = d*(f*dsj+s*dfj)
      else
        s   = a*s*f
        dsi = a*(f*dsi+s*dfi)
        dsj = a*(f*dsj+s*dfj)
      endif
#endif
c
      end subroutine wfn1_entropy2_ab
c
c-----------------------------------------------------------------------
c
!> \brief Entropy enhancement toward low occupations
!>
!> From Berylium atom we know that low occupied orbitals need to resist
!> the tendency towards half occupation, otherwise all the occupation
!> of the 2p orbitals will collapse into a single 2p orbital, generating
!> an unphysical state.
!>
!> Here we will enhance the entropy if the maximum occupation of two
!> low occupied orbitals is pushed further down.
!>
      subroutine wfn1_entropy2_ll(pi,pj,f,dfi,dfj)
      implicit none
!>    The occupation number \f$p_i\f$
      double precision, intent(in)  :: pi
!>    The occupation number \f$p_j\f$
      double precision, intent(in)  :: pj
!>    The enhancement factor \f$f\f$
      double precision, intent(out) :: f
!>    The gradient \f$\partial f/\partial p_i\f$
      double precision, intent(out) :: dfi
!>    The gradient \f$\partial f/\partial p_j\f$
      double precision, intent(out) :: dfj
c
c     Local
c
      double precision :: r   ! the "maximum" of pi and pj
      double precision :: dri ! the gradient of r wrt pi
      double precision :: drj ! the gradient of r wrt pj
      double precision, parameter :: fac0 = 1.0d0
      double precision, parameter :: pow0 = 2.0d0
      double precision, parameter :: fac1 = 1.0d1
      double precision, parameter :: pow1 = 2.0d0*pow0
      double precision, parameter :: fac2 = 1.0d2
      double precision, parameter :: pow2 = 2.0d0*pow1
      double precision :: x
      double precision :: dfx
      double precision :: dxr
c
c     Code
c
cDEBUG
      f   = 1.0d0
      dfi = 0.0d0
      dfj = 0.0d0
      return
cDEBUG
      call wfn1_cont_max(pi,pj,6.0d0,r,dri,drj)
      x   =  0.5d0-r
      dxr = -1.0d0
      f   =  1.0d0
     &    +  fac0*x**pow0
     &    +  fac1*x**pow1
     &    +  fac2*x**pow2
      dfx =  fac0*pow0*x**(pow0-1.0d0)
     &    +  fac1*pow1*x**(pow1-1.0d0)
     &    +  fac2*pow2*x**(pow2-1.0d0)
      dfi =  dfx*dxr*dri
      dfj =  dfx*dxr*drj
      end subroutine wfn1_entropy2_ll
c
c-----------------------------------------------------------------------
c
!> \brief Continuous and continuously differential approximation to max
!>
!> We need a function \f$f(x,y)\f$ that generates a value that is
!> close to \f$max(x,y)\f$ but the function should be continuous and
!> continuously differentiable. Also as \f$x\f$ and \f$y\f$ are orbital
!> occupation numbers their values are limited according to
!> \f$0 \le x,y \le 1\f$.
!>
!> To implement such a function we make use of the fact that the p-norm
!> of a vector for \f$p > 1\f$ is dominated by the largest number in 
!> the vector. So we approximate the maximum function as
!> \f{eqnarray}{
!>   max(x,y) \approx \sqrt[p]{\frac{x^p+y^p}{2}}
!> \f}
!>
      subroutine wfn1_cont_max(x,y,p,r,drx,dry)
      implicit none
!>    The \f$x\f$-argument
      double precision, intent(in)  :: x
!>    The \f$y\f$-argument
      double precision, intent(in)  :: y
!>    The power \f$p\f$ defining the order of the norm
      double precision, intent(in)  :: p
!>    The result \f$r\f$
      double precision, intent(out) :: r
!>    The derivative \f$\partial r/\partial x\f$
      double precision, intent(out) :: drx
!>    The derivative \f$\partial r/\partial y\f$
      double precision, intent(out) :: dry
c
      r   = ((x**p+y**p)/2.0d0)**(1.0d0/p)
      drx = (((x**p+y**p)/2.0d0)**(1.0d0/p-1.0d0))
     &      *(x**(p-1.0d0))/2.0d0
      dry = (((x**p+y**p)/2.0d0)**(1.0d0/p-1.0d0))
     &      *(y**(p-1.0d0))/2.0d0
      end subroutine wfn1_cont_max
c
c-----------------------------------------------------------------------
c
!> \brief Continuous and continuously differential approximation to min
!>
!> We need a function \f$f(x,y)\f$ that generates a value that is
!> close to \f$min(x,y)\f$ but the function should be continuous and
!> continuously differentiable. Also as \f$x\f$ and \f$y\f$ are orbital
!> occupation numbers their values are limited according to
!> \f$0 \le x,y \le 1\f$.
!>
!> To implement such a function we make use of the fact that the p-norm
!> of a vector for \f$p > 1\f$ is dominated by the largest number in 
!> the vector. So we approximate the minimum function as
!> \f{eqnarray}{
!>   min(x,y) \approx \sqrt[-p]{\frac{x^{-p}+y^{-p}}{2}}
!> \f}
!>
      subroutine wfn1_cont_min(x,y,p,r,drx,dry)
      implicit none
!>    The \f$x\f$-argument
      double precision, intent(in)  :: x
!>    The \f$y\f$-argument
      double precision, intent(in)  :: y
!>    The power \f$p\f$ defining the order of the norm
      double precision, intent(in)  :: p
!>    The result \f$r\f$
      double precision, intent(out) :: r
!>    The derivative \f$\partial r/\partial x\f$
      double precision, intent(out) :: drx
!>    The derivative \f$\partial r/\partial y\f$
      double precision, intent(out) :: dry
c
      double precision :: q
c
      q   = -p
      r   = ((x**q+y**q)/2.0d0)**(1.0d0/q)
      drx = (((x**q+y**q)/2.0d0)**(1.0d0/q-1.0d0))
     &      *(x**(q-1.0d0))/2.0d0
      dry = (((x**q+y**q)/2.0d0)**(1.0d0/q-1.0d0))
     &      *(y**(q-1.0d0))/2.0d0
      end subroutine wfn1_cont_min
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_entropy2_aa(pi,pj,pij,s,dsi,dsj,dsij)
      use wfn1_fac
      implicit none
C>    The occupation number of orbital \f$i\f$
      double precision, intent(in)  :: pi
C>    The occupation number of orbital \f$j\f$
      double precision, intent(in)  :: pj
C>    The occupation number of orbital pair \f$i\f$, \f$j\f$
      double precision, intent(in)  :: pij
C>    The entropy
      double precision, intent(out) :: s
C>    The derivative wrt \f$p_i\f$
      double precision, intent(out) :: dsi
C>    The derivative wrt \f$p_j\f$
      double precision, intent(out) :: dsj
C>    The derivative wrt \f$p_{ij}\f$
      double precision, intent(out) :: dsij
c
c     Local
c
      double precision :: a, ai, aj, aij
      double precision ::    si, sj, sij
      double precision :: ds1_i
      double precision :: ds1_j
      double precision :: ds1_ij
      double precision :: ds2_i
      double precision :: ds2_j
c
c     Global
c
c     double precision :: fac_aa
c     double precision :: fac_ab
c     double precision :: trm_aa
c     double precision :: trm_ab
c     double precision :: pow_y_aa
c     double precision :: pow_y_ab
c     double precision :: pow_z_aa
c     double precision :: pow_z_ab
c     double precision, allocatable :: orb_ea(:)
c     double precision, allocatable :: orb_eb(:)
c     common/wfn1_fac/fac_aa,fac_ab,
c    +                trm_aa,trm_ab,
c    +                pow_y_aa,pow_y_ab,
c    +                pow_z_aa,pow_z_ab,
c    +                orb_ea,orb_eb
c
#ifdef OLD
#ifdef LOGGAMMA
      call wfn1_entropy_si(0.5d0,ai,ds1_i)
      call wfn1_entropy_si(0.5d0,aj,ds1_j)
      call wfn1_entropy_sij(0.5d0,0.5d0,aij,ds2_i,ds2_j)
#else
c     Computing the scale factor A from this ensures that the
c     aa and ab entropies have the same scale.
      call wfn1_entropy_si(pow_aa,0.5d0,ai,ds1_i)
      call wfn1_entropy_si(pow_aa,0.5d0,aj,ds1_j)
      call wfn1_entropy_si(pow_aa,0.5d0*0.5d0,aij,ds2_i)
#endif
      a = 0.25d0/(ai+aj-aij)
c
c     The following factor comes from a H2 calculation in the STO-3G
c     basis set at an interatomic distance of the equilibrium bond-
c     length plus 5 Angstrom, and then dividing the UHF 2-electron
c     energy by the RHF 2-electron energy and subtracting that fraction
c     (0.21364) from 1. I.e. electron correlation does not annihilate
c     the 2-electron interaction completely but mostly.
c
c     a = a*0.78636d0*0.83606d0
c     a = a*0.78636d0
c     if (pi.lt.0.5d0.and.pj.lt.0.5d0) then
c       a = fac_aa_ll*a
c     else if (pi.lt.0.5d0.and.pj.ge.0.5d0) then
c       a = fac_aa_hl*a
c     else if (pi.ge.0.5d0.and.pj.lt.0.5d0) then
c       a = fac_aa_hl*a
c     else
c       a = fac_aa_hh*a
c     endif
      a = fac_aa*a
c
      call wfn1_entropy_si(pow_aa,pi,si,ds1_i)
      call wfn1_entropy_si(pow_aa,pj,sj,ds1_j)
      call wfn1_entropy_si(pow_aa,pij,sij,ds1_ij)
c
      s   = a*(si+sj-sij)
cDEBUG
c     if (s.gt.1.0d0) then
c       write(*,*)'S too large: ',s,pi,pj
c     else if (s.lt.0.0d0) then
c       write(*,*)'S too small: ',s,pi,pj
c     endif
cDEBUG
      dsi  =  a*ds1_i
      dsj  =  a*ds1_j
      dsij = -a*ds1_ij
#else
      s    = 0.0d0
      dsi  = 0.0d0
      dsj  = 0.0d0
      dsij = 0.0d0
#endif
c
      end subroutine wfn1_entropy2_aa
c
c-----------------------------------------------------------------------
c
c Note: x*ln(x) is the Sterling approximation to ln(Gamma(x+1)) which 
c       the generalization of ln(x!) to real numbers. The latter is
c       referred to as the log-gamma function and its derivative is
c       called the digamma function. The Sterling approximation is very
c       accurate for large x, but we are interested in 0 <= x <= 1.
c       In that domain the Sterling approximation is not very accurate.
c       Worse the digamma function has a finite value at x=0 whereas the
c       derivative of the Stirling approximation has an asymptote at
c       x=0. The latter makes optimizations very hard. Nevertheless
c       we cannot carelessly switch to the log-gamma and digamma 
c       functions as there are few implementations of these functions.
c       The most commonly available one resides in the GNU scientific
c       library (GSL) which is released under the GPL. Hence we cannot
c       use that implementation in NWChem. So we are between a rock and
c       a hard place here.
c
c-----------------------------------------------------------------------
