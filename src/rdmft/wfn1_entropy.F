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
      subroutine wfn1_entropy1(p,s,ds)
      implicit none
!>    The occupation number
      double precision, intent(in)  :: p
!>    The entropy
      double precision, intent(out) :: s
!>    The entropy gradient
      double precision, intent(out) :: ds
c
      double precision :: a
      double precision :: tol
c
#ifdef USE_GSL
      double precision, external :: gsl_sf_lngamma
      double precision, external :: gsl_sf_psi
      a = log(0.5d0)/(2.0d0*gsl_sf_lngamma(1.5d0))
      s = a*(gsl_sf_lngamma(1.0d0+p) + gsl_sf_lngamma(2.0d0-p))
      ds = a*(gsl_sf_psi(1.0d0+p) - gsl_sf_psi(2.0d0-p))
#else
      tol = 1.0d-10
      if   (p.lt.tol) then
        s  = 0.0d0
        ds = log(tol)
      else (p.gt.1.0d0-tol) then
        s  = 0.0d0
        ds = -log(tol)
      else
        s  = p*log(p) + (1.0d0-p)*log(1.0d0-p)
        ds = log(p) - log(1.0d0-p)
      endif
#endif
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine wfn1_entropy2(p1,p2,s,ds1,ds2)
      implicit none
C>
C> \brief Calculate the 2-electron entropy associated with an 
C> orbital pair and its derivatives
C>
C> The entropy of a natural orbital pair with occupation numbers
C> \f$p_1\f$ and \f$p_2\f$ is
C> given by
C> \f{eqnarray}{
C>    \f$S^2 &=& p_1 p_2\ln(p_1 p_2) 
C>               +(1-p_1)p_2\ln((1-p_1)p_2) \\\\
C>           &&  +p_1(1-p_2)\ln(p_1(1-p_2))
C>               +(1-p_1)(1-p_2)\ln((1-p_1)(1-p_2))
C> \f}.
C> 
      double precision, intent(in)  :: p1  ! The occupation number of
                                           ! orbital 1.
      double precision, intent(in)  :: p2  ! The occupation number of
                                           ! orbital 2.
      double precision, intent(out) :: s   ! The entropy
      double precision, intent(out) :: ds1 ! The entropy derivative for
                                           ! orbital 1.
      double precision, intent(out) :: ds2 ! The entropy derivative for
                                           ! orbital 2.
c
      s =   p1*p2*log(p1*p2)
     &  +   (1.0d0-p1)*p2*log((1.0d0-p1)*p2)
     &  +   p1*(1.0d0-p2)*log(p1*(1.0d0-p2))
     &  +   (1.0d0-p1)*(1.0d0-p2)*log((1.0d0-p1)*(1.0d0-p2))
      ds1 = p2*log(p1*p2) + (2*p2-1.0d0)
     &    - (1.0d0-p2)*log((1.0d0-p1)*(1.0d0-p2))
      ds2 = p1*log(p1*p2) + (2*p1-1.0d0)
     &    - (1.0d0-p1)*log((1.0d0-p1)*(1.0d0-p2))
c
      end
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
