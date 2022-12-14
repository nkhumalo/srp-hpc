#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if !defined(NWAD_PRINT)
C> \ingroup nwxc
C> @{
C>
C> \file nwxc_c_m11.F
C> Implementation of the M08 and M11 correlation functionals
C>
C> @}
#endif
#endif
C>
C> \ingroup nwxc_priv
C> @{
C>
C> \brief The M08 and M11 correlation functionals
C>
C> This routine implements correlation meta-GGA functionals 
C> of the M08 [1] and M11 [2,3] families of functionals.
C>
C> ### References ###
C>
C> [1] Y. Zhao, D. G. Truhlar,
C>     "Exploring the limit of accuracy of the global hybrid meta 
C>     density functional for main-group thermochemistry, kinetics,
C>     and noncovalent interactions",
C>     J. Chem. Theory Comput. <b>4</b>, 1849-1868 (2008), DOI:
C>     <a href="https://doi.org/10.1021/ct800246v">
C>     10.1021/ct800246v</a>.
C>
C> [2] R. Peverati, D. G. Truhlar,
C>     "Improving the accuracy of hybrid meta-GGA density functionals
C>     by range separation",
C>     J. Phys. Chem. Lett. <b>2</b>, 2810-2817 (2011), DOI:
C>     <a href="https://doi.org/10.1021/jz201170d">
C>     10.1021/jz201170d</a>.
C>
C> [3] R. Peverati, D. G. Truhlar,
C>     "M11-L: A local density functional that provides improved
C>     accuracy for electronic structure calculations in chemistry
C>     and physics",
C>     J. Phys. Chem. Lett. <b>3</b>, 117-124 (2011), DOI:
C>     <a href="https://doi.org/10.1021/jz201525m">
C>     10.1021/jz201525m</a>.
C>
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
      Subroutine nwxc_c_m11_p(param, tol_rho, ipol, nq, wght, rho,
     &                        rgamma, tau, func)
#else
      Subroutine nwxc_c_m11(param, tol_rho, ipol, nq, wght, rho, rgamma,
     &                      tau, func)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      Subroutine nwxc_c_m11_d2(param, tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func)
#else
      Subroutine nwxc_c_m11_d3(param, tol_rho, ipol, nq, wght, 
     &                         rho, rgamma, tau, func)
#endif

c
c$Id$
c
c
c**********************************************************************c
c                                                                      c
c  xc_cm11 evaluates the correlation part of the M08 and M11 suite of  c
c  functionals on the grid.                                            c
c  !!! Second derivatives are not available yet.                       c
c                                                                      c
c  Ref: (a) Zhao, Y.  and Truhlar, D. G. JCTC, 2008, 4 , 1849          c
c       (b) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 2, 2810  c
c       (c) Peverati, R. and Truhlar, D. G. J.P.C.Lett. 2011, 3, 117   c
c                                                                      c
c       ijzy - 1 M08-HX (a)                                            c
c       ijzy - 2 M08-SO (a)                                            c
c       ijzy - 3 M11 (b)                                               c
c       ijzy - 4 M11-L (c)                                             c
c                                                                      c
c Coded by Roberto Peverati (12/11)                                    c
c                                                                      c
c**********************************************************************c
c
#include "nwad.fh"
c
      implicit none
c
#include "nwxc_param.fh"
c
c     Input and other parameters
c
#if defined(NWAD_PRINT)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
      type(nwad_dble)::param(*)
      type(nwad_dble)::at0, at1, at2, at3, at4, at5
      type(nwad_dble)::at6, at7, at8, at9, at10, at11
      type(nwad_dble)::bt0, bt1,  bt2, bt3, bt4, bt5
      type(nwad_dble)::bt6, bt7, bt8, bt9, bt10, bt11
#else
      double precision param(*)
      double precision at0, at1, at2, at3, at4, at5
      double precision at6, at7, at8, at9, at10, at11
      double precision bt0, bt1,  bt2, bt3, bt4, bt5
      double precision bt6, bt7, bt8, bt9, bt10, bt11
#endif
#else
      double precision param(*) !< [Input] Parameters of the functional
                                !< (see Table 2 of [1], and Table 1 of [2])
                                !< - param(1): \f$ a_0 \f$
                                !< - param(2): \f$ a_1 \f$
                                !< - param(3): \f$ a_2 \f$
                                !< - param(4): \f$ a_3 \f$
                                !< - param(5): \f$ a_4 \f$
                                !< - param(6): \f$ a_5 \f$
                                !< - param(7): \f$ a_6 \f$
                                !< - param(8): \f$ a_7 \f$
                                !< - param(9): \f$ a_8 \f$
                                !< - param(10): \f$ a_9 \f$
                                !< - param(11): \f$ a_{10} \f$
                                !< - param(12): \f$ a_{11} \f$
                                !< - param(13): \f$ b_0 \f$
                                !< - param(14): \f$ b_1 \f$
                                !< - param(15): \f$ b_2 \f$
                                !< - param(16): \f$ b_3 \f$
                                !< - param(17): \f$ b_4 \f$
                                !< - param(18): \f$ b_5 \f$
                                !< - param(19): \f$ b_6 \f$
                                !< - param(20): \f$ b_7 \f$
                                !< - param(21): \f$ b_8 \f$
                                !< - param(22): \f$ b_9 \f$
                                !< - param(23): \f$ b_{10} \f$
                                !< - param(24): \f$ b_{11} \f$
      double precision at0, at1, at2, at3, at4, at5
      double precision at6, at7, at8, at9, at10, at11
      double precision bt0, bt1,  bt2, bt3, bt4, bt5
      double precision bt6, bt7, bt8, bt9, bt10, bt11
#endif
      double precision tol_rho !< [Input] The lower limit on the density
      integer nq               !< [Input] The number of points
      integer ipol             !< [Input] The number of spin channels
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density
c
      type(nwad_dble)::rho(nq,*) !< [Input] The density
c
c     Charge Density Gradient Norm
c
      type(nwad_dble)::rgamma(nq,*) !< [Input] The density gradient norm
c
c     Kinetic Energy Density
c
      type(nwad_dble)::tau(nq,*) !< [Input] The kinetic energy density
c
c     Functional values
c
      type(nwad_dble)::func(*) !< [Output] The functional value
c
c     Sampling Matrices for the XC Potential
c
c     double precision Amat(nq,*) !< [Output] Derivative wrt density
c     double precision Cmat(nq,*) !< [Output] Derivative wrt rgamma
c     double precision Mmat(nq,*) !< [Output] Derivative wrt tau
c
c     Constants in PBE functional
c
      double precision GAMMA, BETA, PI
      parameter (GAMMA = 0.03109069086965489503494086371273d0)
      parameter (BETA = 0.06672455060314922d0)
      parameter (PI = 3.1415926535897932385d0)
c
c     Threshold parameters
c
      double precision EXPTOL
      parameter (EXPTOL = 40.0d0)
      type(nwad_dble)::rho_t(3)
c     type(nwad_dble)::dsqgamma
c
c     Intermediate derivative results, etc.
c
      integer n
      
      type(nwad_dble)::rhoa, rhob, rhoval, gammaval, TA, TB
      type(nwad_dble)::E1, E2
      double precision de1drb, de2dga, de2dgb, de2dgc, de2dk, de2dra
      double precision de2drb, df1dr, df1dtau, df1dw, df2dr, df2dtau
      double precision de2dgaa, de2dgab, de2dgbb
      double precision de1dk, de1dra
      double precision df2dw, dggadga, dggadgb, dggadgc, dggadra
      double precision dggadrb, dgrhody, dhdga, dhdgb, dhdgc, dhdy
      double precision dlsdadra, dlsdadrb, dtdr, dtdtau, dwdt
      type(nwad_dble)::fsig1, fsig2
c     type(nwad_dble)::taup, tauueg, tsig, wsig 
      type(nwad_dble)::taup, tauueg, wsig 
      
      type(nwad_dble)::nepsc(1),epsc
      double precision dnepscdn(2)
      double precision depscdna, depscdnb
      type(nwad_dble)::H0
      double precision dH0dna, dH0dnb, dH0dg

      type(nwad_dble)::phi,zeta,arglog,fAt,fAtnum,fAtden,A,t2,ks
      type(nwad_dble)::argexp
      double precision dphidna, dphidnb, dphidzeta
      double precision dzetadna, dzetadnb
      double precision darglogdna, darglogdnb, darglogdg
      double precision dfAtdt, dfAtdA
      double precision dfAtnumdt, dfAtnumdA
      double precision dfAtdendt, dfAtdendA
      double precision dfAtdna, dfAtdnb, dfAtdg
      double precision dAdna, dAdnb
      double precision dtdna, dtdnb, dtdg
      double precision dksdna, dksdnb
      double precision dargexpdna, dargexpdnb
      type(nwad_dble)::expinA
    
      type(nwad_dble)::LSDA,GGA
 
      double precision f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
      double precision f1o3, f2o3, f4o3, f5o3     
      Save F1, F2, F3, F4, F5, F6, F7, F8, F9, F10, F11
      DATA F1/1.0D+00/,  F2/2.0D+00/,  F3/3.0D+00/,
     $     F4/4.0D+00/,  F5/5.0D+00/,  F6/6.0D+00/,
     $     F7/7.0D+00/,  F8/8.0D+00/,  F9/9.0D+00/, 
     $     F10/10.0D+00/,F11/11.0D+00/
c
       at0=param( 1)
       at1=param( 2)
       at2=param( 3)
       at3=param( 4)
       at4=param( 5)
       at5=param( 6)
       at6=param( 7)
       at7=param( 8)
       at8=param( 9)
       at9=param(10)
       at10=param(11)
       at11=param(12)
c
       bt0=param(13)
       bt1=param(14)
       bt2=param(15)
       bt3=param(16)
       bt4=param(17)
       bt5=param(18)
       bt6=param(19)
       bt7=param(20)
       bt8=param(21)
       bt9=param(22)
       bt10=param(23)
       bt11=param(24)
c
c       at0=0.000000D+00
c       at1=0.000000D+00
c       at2=0.000000D+00
c       at3=0.000000D+00
c       at4=0.000000D+00
c       at5=0.000000D+00
c       at6=0.000000D+00
c       at7=0.000000D+00
c       at8=0.000000D+00
c       at9=0.000000D+00
c       at10=0.000000D+00 
c       at11=0.000000D+00 
c       bt0=0.000000D+00
c       bt1=0.000000D+00
c       bt2=0.000000D+00
c       bt3=0.000000D+00
c       bt4=0.000000D+00
c       bt5=0.000000D+00
c       bt6=0.000000D+00
c       bt7=0.000000D+00
c       bt8=0.000000D+00
c       bt9=0.000000D+00
c       bt10=0.000000D+00 
c       bt11=0.000000D+00 
c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c

c     if (ijzy.eq.1) then
c     Parameters for M08-HX
c       at0=    1.0000000D+00
c       at1=   -4.0661387D-01
c       at2=   -3.3232530D+00
c       at3=    1.5540980D+00
c       at4=    4.4248033D+01
c       at5=   -8.4351930D+01
c       at6=   -1.1955581D+02
c       at7=    3.9147081D+02
c       at8=    1.8363851D+02
c       at9=   -6.3268223D+02
c       at10=  -1.1297403D+02
c       at11=   3.3629312D+02

c       bt0=    1.3812334D+00
c       bt1=   -2.4683806D+00
c       bt2=   -1.1901501D+01
c       bt3=   -5.4112667D+01
c       bt4=    1.0055846D+01
c       bt5=    1.4800687D+02
c       bt6=    1.1561420D+02
c       bt7=    2.5591815D+02
c       bt8=    2.1320772D+02
c       bt9=   -4.8412067D+02
c       bt10=  -4.3430813D+02
c       bt11=   5.6627964D+01
c      elseif (ijzy.eq.2) then
c     Parameters for M08-SO 
c       at0=    1.0000000D+00
c       at1=    0.0000000D+00
c       at2=   -3.9980886D+00
c       at3=    1.2982340D+01
c       at4=    1.0117507D+02
c       at5=   -8.9541984D+01
c       at6=   -3.5640242D+02
c       at7=    2.0698803D+02
c       at8=    4.6037780D+02
c       at9=   -2.4510559D+02
c       at10=  -1.9638425D+02
c       at11=   1.1881459D+02

c       bt0=    1.0000000D+00
c       bt1=   -4.4117403D+00
c       bt2=   -6.4128622D+00
c       bt3=    4.7583635D+01
c       bt4=    1.8630053D+02
c       bt5=   -1.2800784D+02
c       bt6=   -5.5385258D+02
c       bt7=    1.3873727D+02
c       bt8=    4.1646537D+02
c       bt9=   -2.6626577D+02
c       bt10=   5.6676300D+01
c       bt11=   3.1673746D+02
c      elseif (ijzy.eq.3) then
c     Parameters for M11
c       at0=   1.0000000D+00
c       at1=   0.0000000D+00
c       at2=  -3.8933250D+00
c       at3=  -2.1688455D+00
c       at4=   9.3497200D+00
c       at5=  -1.9845140D+01
c       at6=   2.3455253D+00
c       at7=   7.9246513D+01
c       at8=   9.6042757D+00
c       at9=  -6.7856719D+01
c       at10= -9.1841067D+00
c       at11=  0.0000000D+00

c       bt0=   7.2239798D-01
c       bt1=   4.3730564D-01
c       bt2=  -1.6088809D+01
c       bt3=  -6.5542437D+01
c       bt4=   3.2057230D+01
c       bt5=   1.8617888D+02
c       bt6=   2.0483468D+01
c       bt7=  -7.0853739D+01
c       bt8=   4.4483915D+01
c       bt9=  -9.4484747D+01
c       bt10= -1.1459868D+02
c       bt11=  0.0000000D+00
c      elseif (ijzy.eq.4) then
c     Parameters for M11-L
c       at0=   1.000000D+00
c       at1=   0.000000D+00
c       at2=   2.750880D+00
c       at3=  -1.562287D+01
c       at4=   9.363381D+00
c       at5=   2.141024D+01
c       at6=  -1.424975D+01
c       at7=  -1.134712D+01
c       at8=   1.022365D+01
c       at9=   0.000000D+00
c       at10=  0.000000D+00
c       at11=  0.000000D+00
c
c       bt0=   1.000000D+00
c       bt1=  -9.082060D+00
c       bt2=   6.134682D+00
c       bt3=  -1.333216D+01
c       bt4=  -1.464115D+01
c       bt5=   1.713143D+01
c       bt6=   2.480738D+00
c       bt7=  -1.007036D+01
c       bt8=  -1.117521D-01
c       bt9=   0.000000D+00
c       bt10=  0.000000D+00
c       bt11=  0.000000D+00
c      endif


      F1o3 = F1/F3
      F2o3 = F2/F3
      F4o3 = F4/F3
      F5o3 = F5/F3
      


c
c     ======> BOTH SPIN-RESTRICETED AND UNRESTRICTED <======
c
      do 20 n = 1, nq

      if (ipol.eq.1) then
        RHOA = rho(n,R_T)/F2
        RHOB = RHOA
        rho_t(1) = rho(n,R_T)        
c       gammaval =(delrho(n,1,1)*delrho(n,1,1) +
c    &             delrho(n,2,1)*delrho(n,2,1) +
c    &             delrho(n,3,1)*delrho(n,3,1))
        gammaval = rgamma(n,G_TT)
        TA = tau(n,T_T)/F2
        TB = TA
      else
        RHOA = rho(n,R_A)
        RHOB = rho(n,R_B)
c       gammaval = delrho(n,1,1)*delrho(n,1,1) +
c    &             delrho(n,1,2)*delrho(n,1,2) +
c    &             delrho(n,2,1)*delrho(n,2,1) +
c    &             delrho(n,2,2)*delrho(n,2,2) +
c    &             delrho(n,3,1)*delrho(n,3,1) +
c    &             delrho(n,3,2)*delrho(n,3,2) +
c    &       2.d0*(delrho(n,1,1)*delrho(n,1,2) +
c    &             delrho(n,2,1)*delrho(n,2,2) +
c    &             delrho(n,3,1)*delrho(n,3,2))
        if (rhoa.gt.0.5d0*tol_rho.and.rhob.gt.0.5d0*tol_rho) then
          gammaval = rgamma(n,G_AA)+rgamma(n,G_BB)+2.0d0*rgamma(n,G_AB)
        elseif (rhoa.gt.0.5d0*tol_rho.and.rhob.le.0.5d0*tol_rho) then
          rhob = 0.0d0
          gammaval = rgamma(n,G_AA)
        elseif (rhoa.le.0.5d0*tol_rho.and.rhob.gt.0.5d0*tol_rho) then
          rhoa = 0.0d0
          gammaval = rgamma(n,G_BB)
        endif
        rho_t(1) = rhoa
        rho_t(2) = rhob
        TA =  tau(n,T_A)
        TB =  tau(n,T_B)
      endif
      RHOVAL  = RHOA + RHOB
      TauP = TA+TB
c     t = tol_rho
c     dsqgamma = max(sqrt(gammaval),t)
c     dsqgamma = sqrt(gammaval)
      if (rhoval.le.tol_rho) goto 20
c     if (taup.lt.tol_rho) goto 20 
c
      nepsc(1) = 0.0d0
      dnepscdn(1) = 0.0d0
      if (ipol.eq.2) dnepscdn(2) = 0.0d0
c
c        call for LDA bit
c        this implementation temporarily assigns the pw91LDA for
c        use in the metaGGA local part
c
c           call xc_pw91lda(tol_rho,1d0,.true.,.false.,rho_t,
c    &         dnepscdn,1,ipol,nepsc,1d0,
c    &         .false.,dummy)
#if !defined(SECOND_DERIV) && !defined(THIRD_DERIV)
#if defined(NWAD_PRINT)
         call nwxc_c_pw91lda_p(tol_rho,ipol,1,1.0d0,rho_t,
     &        nepsc)
#else
         call nwxc_c_pw91lda(tol_rho,ipol,1,1.0d0,rho_t,
     &        nepsc)
#endif
#elif defined(SECOND_DERIV) && !defined(THIRD_DERIV)
         call nwxc_c_pw91lda_d2(tol_rho,ipol,1,1.0d0,rho_t,
     &        nepsc)
#else
         call nwxc_c_pw91lda_d3(tol_rho,ipol,1,1.0d0,rho_t,
     &        nepsc)
#endif
        
c
c        PBE non-local part
c        This part is a modified version of xc_cMpbe96.F
c
         if(abs(nepsc(1)).lt.tol_rho*tol_rho) goto 20
c
         epsc = nepsc(1)/rhoval
c        if (ipol.eq.1) then
c           depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
c           depscdnb = depscdna
c        else
c           depscdna = dnepscdn(1)/rhoval-nepsc/(rhoval**2)
c           depscdnb = dnepscdn(2)/rhoval-nepsc/(rhoval**2)
c        endif
c
         ks = 2.0d0*(3.0d0*PI*PI*rhoval)**(1.0d0/6.0d0)/dsqrt(PI)
c        dksdna = (1.0d0/6.0d0)*ks/rhoval
c        dksdnb = dksdna
c
         if (ipol.eq.1) then
            zeta = 0.0d0
         else
            zeta = (rho(n,R_A)-rho(n,R_B))/rhoval
         endif
         if(zeta.lt.-1.0d0) zeta=-1.0d0
         if(zeta.gt. 1.0d0) zeta= 1.0d0
         phi = 0.0d0
         if (1.0d0+zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*(1.0d0+zeta)**(2.0d0/3.0d0)
         endif
         if (1.0d0-zeta.gt.1.0d-10) then
           phi = phi + 0.5d0*(1.0d0-zeta)**(2.0d0/3.0d0)
         endif
c        if (ipol.eq.1) then
c           dzetadna = 1.0d0/rhoval
c           dzetadnb = -1.0d0/rhoval
c        else
c           dzetadna =  2.0d0*rho(n,R_B)/(rhoval**2)
c           dzetadnb = -2.0d0*rho(n,R_A)/(rhoval**2)
c        endif
c        phi = 0.5d0*((1.0d0+zeta)**(2.0d0/3.0d0)
c    &               +(1.0d0-zeta)**(2.0d0/3.0d0))
c        if (dabs(1.0d0-zeta).lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &             (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta))
c        else if (dabs(1.0d0+zeta).lt.tol_rho) then
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &            -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
c        else
c           dphidzeta = 0.5d0*(2.0d0/3.0d0)*(
c    &         (1.0d0+zeta)**(2.0d0/3.0d0)/(1.0d0+zeta)
c    &        -(1.0d0-zeta)**(2.0d0/3.0d0)/(1.0d0-zeta))
c        endif
c        dphidna = dphidzeta*dzetadna
c        dphidnb = dphidzeta*dzetadnb
c
c        t = dsqgamma/(2.0d0*phi*ks*rhoval)
         t2 = gammaval/((2.0d0*phi*ks*rhoval)**2)
c        dtdna = -t/rhoval-t/phi*dphidna-t/ks*dksdna
c        dtdnb = -t/rhoval-t/phi*dphidnb-t/ks*dksdnb
c
         argexp = -epsc/GAMMA/(phi**3)
c        dargexpdna = -depscdna/GAMMA/(phi**3)
c    &                +3.0d0*epsc/GAMMA/(phi**4)*dphidna
c        dargexpdnb = -depscdnb/GAMMA/(phi**3)
c    &                +3.0d0*epsc/GAMMA/(phi**4)*dphidnb
c
         if (abs(argexp).lt.EXPTOL) then
            expinA=exp(argexp)
         else
            expinA=0.0d0
         endif
         A = BETA/GAMMA/(expinA-1.0d0)
c        dAdna = -BETA/GAMMA*dargexpdna*expinA/(expinA-1.0d0)**2
c        dAdnb = -BETA/GAMMA*dargexpdnb*expinA/(expinA-1.0d0)**2
c
c        fAtnum = 1.0d0+A*t**2.0d0
c        fAtden = 1.0d0+A*t**2.0d0+A**2.0d0*t**4.0d0
         fAtnum = 1.0d0+A*t2
         fAtden = 1.0d0+A*t2+A**2*t2**2
         fAt = fAtnum/fAtden
c        dfAtnumdt = 2.0d0*A*t
c        dfAtnumdA = t**2
c        dfAtdendt = 2.0d0*A*t+4.0d0*A**2*t**3
c        dfAtdendA = t**2+2.0d0*A*t**4
c        dfAtdt = (dfAtnumdt*fAtden-fAtnum*dfAtdendt)/(fAtden**2)
c        dfAtdA = (dfAtnumdA*fAtden-fAtnum*dfAtdendA)/(fAtden**2)
c        dfAtdna = dfAtdt * dtdna + dfAtdA * dAdna
c        dfAtdnb = dfAtdt * dtdnb + dfAtdA * dAdnb
c
c        arglog = 1.0d0 + BETA/GAMMA*t**2.0d0*fAt
         arglog = 1.0d0 + BETA/GAMMA*t2*fAt
c        darglogdna = BETA/GAMMA*(2.0d0*t*dtdna*fAt
c    &                            +t*t*dfAtdna)
c        darglogdnb = BETA/GAMMA*(2.0d0*t*dtdnb*fAt
c    &                            +t*t*dfAtdnb)
c
         H0 = GAMMA*(phi**3)*log(arglog)
c        dH0dna = GAMMA*(3.0d0*(phi**2)*dphidna*dlog(arglog)
c    &                  +(phi**3)*darglogdna/arglog)
c        dH0dnb = GAMMA*(3.0d0*(phi**2)*dphidnb*dlog(arglog)
c    &                  +(phi**3)*darglogdnb/arglog)       
c
        TauUEG=F3*(F3*Pi*Pi)**F2o3*rhoval**F5o3/F10
c       Tsig =TauUEG/TauP
c       Wsig =(Tsig - F1)/(Tsig + F1)
        Wsig =(TauUEG - TauP)/(TauUEG + TauP)
        Fsig1=(at0 + Wsig*(at1 + Wsig*(at2 + Wsig*(at3 + Wsig*(
     &            at4 + Wsig*(at5 + Wsig*(at6 + Wsig*(at7 + Wsig*(
     &            at8 + Wsig*(at9 + Wsig*(at10+Wsig*at11)))))))))))

        Fsig2=(bt0 + Wsig*(bt1 + Wsig*(bt2 + Wsig*(bt3 + Wsig*(
     &            bt4 + Wsig*(bt5 + Wsig*(bt6 + Wsig*(bt7 + Wsig*(
     &            bt8 + Wsig*(bt9 + Wsig*(bt10+Wsig*bt11)))))))))))
        LSDA = rhoval*epsc
        GGA  = rhoval*H0
        E1   = LSDA*Fsig1
        E2   = GGA*Fsig2
c
c        derivatives of the meta-GGA part, and total
c
c        dF1dW=( at1 + Wsig*(F2  *at2 + Wsig*(F3*at3 + Wsig*(
c    &            F4 *at4 + Wsig*(F5 *at5 + Wsig*(F6  *at6 + Wsig*(
c    &            F7*at7 + Wsig*(F8*at8 + Wsig*(F9 *at9 + Wsig*(
c    &            F10  *at10+ Wsig*F11*at11))))))))))
c        dF2dW=( bt1 + Wsig*(F2  *bt2 + Wsig*(F3*bt3 + Wsig*(
c    &            F4 *bt4 + Wsig*(F5 *bt5 + Wsig*(F6  *bt6 + Wsig*(
c    &            F7*bt7 + Wsig*(F8*bt8 + Wsig*(F9 *bt9 + Wsig*(
c    &            F10  *bt10+ Wsig*F11*bt11))))))))))
c        dWdT = F2/((F1 + Tsig)**F2)
c        dTdR = Tsig*F5/(F3*rhoval) 
c        dTdTau = -Tsig/TauP
c        dF1dR = dF1dW*dWdT*dTdR
c        dF1dTau=dF1dW*dWdT*dTdTau
c        dF2dR = dF2dW*dWdT*dTdR
c        dF2dTau=dF2dW*dWdT*dTdTau

c        dLSDAdRA=dnepscdn(1)
c        if (ipol.eq.2) dLSDAdRB=dnepscdn(2)
c        dGGAdRA=rhoval*dH0dna+H0
c        if (ipol.eq.2) dGGAdRB=rhoval*dH0dnb+H0
C
c        dE1dRA = dLSDAdRA*Fsig1 + LSDA*dF1dR
c        dE1dRB = 0.0d0
c        if (ipol.eq.2) dE1dRB = dLSDAdRB*Fsig1 + LSDA*dF1dR
C
c        dE2dRA = dGGAdRA*Fsig2 + GGA*dF2dR
c        dE2dRB = 0.0d0
c        if (ipol.eq.2) dE2dRB = dGGAdRB*Fsig2 + GGA*dF2dR 
c
c        dtdg = 0.25d0/(phi*ks*rhoval)/dsqgamma
c        dfAtdg = dfAtdt*dtdg
c        darglogdg = BETA/GAMMA*(2.0d0*t*dtdg*fAt+t*t*dfAtdg)
c        dH0dg = GAMMA*(phi**3)*darglogdg/arglog
C
c        dE2dGAA =    rhoval*dH0dg*Fsig2
c        dE2dGAB = F2*rhoval*dH0dg*Fsig2
c        dE2dGBB =    rhoval*dH0dg*Fsig2
c
c        dE1dK = LSDA*dF1dTau
c        dE2dK = GGA*dF2dTau
c
c        Ec = Ec+ (E1+E2)*qwght(n)
         func(n) = func(n)+(E1+E2)*wght
c        Amat(n,D1_RA) = Amat(n,D1_RA) + (dE1dRA + dE2dRA)*wght
c        if (ipol.eq.2) then
c          Amat(n,D1_RB) = Amat(n,D1_RB) + (dE1dRB+dE2dRB)*wght
c        endif
c        if (ipol.eq.1) then
c          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dE2dGAA*wght
c          Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dE2dGAB*wght
c        else
c          Cmat(n,D1_GAA) = Cmat(n,D1_GAA) + dE2dGAA*wght
c          Cmat(n,D1_GAB) = Cmat(n,D1_GAB) + dE2dGAB*wght
c          Cmat(n,D1_GBB) = Cmat(n,D1_GBB) + dE2dGBB*wght
c        endif
c        Mmat(n,D1_TA)=Mmat(n,D1_TA)+ (de1dk+de2dk)/F2*wght
c        if (ipol.eq.2) then
c          Mmat(n,D1_TB)=Mmat(n,D1_TB)+(de1dk+de2dk)/F2*wght
c        endif
   20 continue
      end
#ifndef NWAD_PRINT
#define NWAD_PRINT
c
c     Compile source again for Maxima
c
#include "nwxc_c_m11.F"
#endif
#ifndef SECOND_DERIV
#define SECOND_DERIV
c
c     Compile source again for the 2nd derivative case
c
#include "nwxc_c_m11.F"
#endif
#ifndef THIRD_DERIV
#define THIRD_DERIV
c
c     Compile source again for the 3rd derivative case
c
#include "nwxc_c_m11.F"
#endif
#undef NWAD_PRINT
C>
C> @}
