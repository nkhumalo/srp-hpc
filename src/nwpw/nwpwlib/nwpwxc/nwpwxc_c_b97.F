C> \ingroup nwpwxc
C> @{
C>
C> \file nwpwxc_c_b97.F
C> The B97 family of correlation functionals
C>
C> @}
#include "nwpwxcP_xc_b97.h"
C>
C> \ingroup nwpwxc_priv
C> @{
C>
C> \brief Evaluate the B97 family of correlation functionals
C>
C> This code evaluates correlation functionals from the 
C> B97 family of functionals [1,2].
C>
C> ### References ###
C>
C> [1] A.D. Becke, "Density-functional thermochemistry. V. Systematic
C>     optimization of exchange-correlation functionals", J. Chem. Phys.
C>     107 (1997) 8554-8560, DOI:
C>     <a href="https://doi.org/10.1063/1.475007">
C>     10.1063/1.475007</a>.
C>
C> [2] S. Grimme, "Semiempirical GGA-type density functional constructed
C>     with a long-range dispersion correction", J. Comput. Chem. 27
C>     (2006) 1787-1799, DOI:
C>     <a href="https://doi.org/10.1002/jcc.20495">
C>     10.1002/jcc.20495</a>.
C>
      Subroutine nwpwxc_c_b97(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                      func,Amat,Cmat)
c
c     $Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision param(*)!< [Input] Parameters of functional as
                               !< defined in [1]:
                               !< - param(1): \f$m\f$ of Eqs.(20).
                               !< - param(2): \f$C_{C\sigma\sigma,0}\f$
                               !< - param(3): \f$C_{C\alpha\beta,0}\f$
                               !< - param(4): \f$C_{C\sigma\sigma,1}\f$
                               !< - param(5): \f$C_{C\alpha\beta,1}\f$
                               !< - param(6): \f$C_{C\sigma\sigma,2}\f$
                               !< - param(7): \f$C_{C\alpha\beta,2}\f$
                               !< - param(8): \f$C_{C\sigma\sigma,3}\f$
                               !< - param(9): \f$C_{C\alpha\beta,3}\f$
                               !< - param(10): \f$C_{C\sigma\sigma,4}\f$
                               !< - param(11): \f$C_{C\alpha\beta,4}\f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)     !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
c
c     Local variables
c
      integer i
      double precision rho_a(0:1)
      double precision rho_b(0:1)
      double precision FC(0:_FXC_NUMDERI)
c
c     Code
c
      if (ipol.eq.1) then
        do i = 1, nq
          rho_a(0) = rho(i,R_T)*0.5d0
          rho_b(0) = rho_a(0)
          rho_a(1) = rgamma(i,G_TT)*0.25d0
          rho_b(1) = rho_a(1)
          if (rho_a(0).gt.tol_rho) then
            call nwpwxcp_c_b97(rho_a,rho_b,wght,tol_rho,FC,param)
            func(i)        = func(i)        + FC(_FXC_E)*wght
            Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
            Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
          endif
        enddo
      else
        do i = 1, nq
          rho_a(0) = rho(i,R_A)
          rho_b(0) = rho(i,R_B)
          rho_a(1) = rgamma(i,G_AA)
          rho_b(1) = rgamma(i,G_BB)
          if (rho_a(0).gt.tol_rho.or.rho_b(0).gt.tol_rho) then
            call nwpwxcp_c_b97(rho_a,rho_b,wght,tol_rho,FC,param)
            func(i)        = func(i)        + FC(_FXC_E)*wght
            Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
            Amat(i,D1_RB)  = Amat(i,D1_RB)  + FC(_FXC_RB)*wght
            Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
            Cmat(i,D1_GBB) = Cmat(i,D1_GBB) + FC(_FXC_GBB)*wght
          endif
        enddo
      endif
c
      end
c
      Subroutine nwpwxc_c_b97_d2(param,tol_rho,ipol,nq,wght,rho,rgamma,
     &                         func,Amat,Amat2,Cmat,Cmat2)
c
c     $Id$
c
      implicit none
c
#include "nwpwxc_param.fh"
c
c     Input and other parameters
c
      double precision param(*)!< [Input] Parameters of functional as
                               !< defined in [1]:
                               !< - param(1): \f$m\f$ of Eqs.(20).
                               !< - param(2): \f$C_{C\sigma\sigma,0}\f$
                               !< - param(3): \f$C_{C\alpha\beta,0}\f$
                               !< - param(4): \f$C_{C\sigma\sigma,1}\f$
                               !< - param(5): \f$C_{C\alpha\beta,1}\f$
                               !< - param(6): \f$C_{C\sigma\sigma,2}\f$
                               !< - param(7): \f$C_{C\alpha\beta,2}\f$
                               !< - param(8): \f$C_{C\sigma\sigma,3}\f$
                               !< - param(9): \f$C_{C\alpha\beta,3}\f$
                               !< - param(10): \f$C_{C\sigma\sigma,4}\f$
                               !< - param(11): \f$C_{C\alpha\beta,4}\f$
      double precision tol_rho !< [Input] The lower limit on the density
      integer ipol             !< [Input] The number of spin channels
      integer nq               !< [Input] The number of points
      double precision wght    !< [Input] The weight of the functional
c
c     Charge Density 
c
      double precision rho(nq,*)    !< [Input] The density
c
c     Charge Density Gradient
c
      double precision rgamma(nq,*) !< [Input] The norm of the density gradients
c
c     Sampling Matrices for the XC Potential
c
      double precision func(nq)     !< [Output] The value of the functional
      double precision Amat(nq,*)   !< [Output] The derivative wrt rho
      double precision Cmat(nq,*)   !< [Output] The derivative wrt rgamma
      double precision Amat2(nq,*)  !< [Output] The 2nd derivative wrt rho
      double precision Cmat2(nq,*)  !< [Output] The 2nd derivative wrt
                                    !< rho and rgamma
c
c     Local variables
c
      integer i
      double precision rho_a(0:1)
      double precision rho_b(0:1)
      double precision FC(0:_FXC_NUMDERI)
c
c     Code
c
      if (ipol.eq.1) then
        do i = 1, nq
          rho_a(0) = rho(i,R_T)*0.5d0
          rho_b(0) = rho_a(0)
          rho_a(1) = rgamma(i,G_TT)*0.25d0
          rho_b(1) = rho_a(1)
          if (rho_a(0).gt.tol_rho) then
            call nwpwxcp_c_b97(rho_a,rho_b,wght,tol_rho,FC,param)
            func(i)        = func(i)        + FC(_FXC_E)*wght
            Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
            Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
            Amat2(i,D2_RA_RA) = Amat2(i,D2_RA_RA) + FC(_FXC_RARA)*wght
            Amat2(i,D2_RA_RB) = Amat2(i,D2_RA_RB) + FC(_FXC_RARB)*wght
            Cmat2(i,D2_RA_GAA) = Cmat2(i,D2_RA_GAA)
     &                         + FC(_FXC_RAGAA)*wght
            Cmat2(i,D2_RA_GBB) = Cmat2(i,D2_RA_GBB)
     &                         + FC(_FXC_RAGBB)*wght
            Cmat2(i,D2_GAA_GAA) = Cmat2(i,D2_GAA_GAA) 
     &                          + FC(_FXC_GAAGAA)*wght
            Cmat2(i,D2_GAA_GBB) = Cmat2(i,D2_GAA_GBB) 
     &                          + FC(_FXC_GAAGBB)*wght
          endif
        enddo
      else
        do i = 1, nq
          rho_a(0) = rho(i,R_A)
          rho_b(0) = rho(i,R_B)
          rho_a(1) = rgamma(i,G_AA)
          rho_b(1) = rgamma(i,G_BB)
          if (rho_a(0).gt.tol_rho.or.rho_b(0).gt.tol_rho) then
            call nwpwxcp_c_b97(rho_a,rho_b,wght,tol_rho,FC,param)
            func(i)        = func(i)        + FC(_FXC_E)*wght
            Amat(i,D1_RA)  = Amat(i,D1_RA)  + FC(_FXC_RA)*wght
            Amat(i,D1_RB)  = Amat(i,D1_RB)  + FC(_FXC_RB)*wght
            Cmat(i,D1_GAA) = Cmat(i,D1_GAA) + FC(_FXC_GAA)*wght
            Cmat(i,D1_GBB) = Cmat(i,D1_GBB) + FC(_FXC_GBB)*wght
            Amat2(i,D2_RA_RA) = Amat2(i,D2_RA_RA) + FC(_FXC_RARA)*wght
            Amat2(i,D2_RA_RB) = Amat2(i,D2_RA_RB) + FC(_FXC_RARB)*wght
            Amat2(i,D2_RB_RB) = Amat2(i,D2_RB_RB) + FC(_FXC_RBRB)*wght
            Cmat2(i,D2_RA_GAA) = Cmat2(i,D2_RA_GAA)
     &                         + FC(_FXC_RAGAA)*wght
            Cmat2(i,D2_RA_GBB) = Cmat2(i,D2_RA_GBB)
     &                         + FC(_FXC_RAGBB)*wght
            Cmat2(i,D2_RB_GAA) = Cmat2(i,D2_RB_GAA)
     &                         + FC(_FXC_RBGAA)*wght
            Cmat2(i,D2_RB_GBB) = Cmat2(i,D2_RB_GBB)
     &                         + FC(_FXC_RBGBB)*wght
            Cmat2(i,D2_GAA_GAA) = Cmat2(i,D2_GAA_GAA) 
     &                          + FC(_FXC_GAAGAA)*wght
            Cmat2(i,D2_GAA_GBB) = Cmat2(i,D2_GAA_GBB) 
     &                          + FC(_FXC_GAAGBB)*wght
            Cmat2(i,D2_GBB_GBB) = Cmat2(i,D2_GBB_GBB) 
     &                          + FC(_FXC_GBBGBB)*wght
          endif
        enddo
      endif
c
      end
c
      subroutine nwpwxcp_c_pwlda(ra,rb,FCLDA)
      implicit none
c
      double precision ra
      double precision rb
      double precision FCLDA(0:_FCLDA_ELEMENTS)
c
      double precision ec
      double precision rho(2)
      double precision Amat(2)
      double precision Amat2(3)
c
      rho(R_A)        = ra
      rho(R_B)        = rb
      ec              = 0.0d0
      Amat(D1_RA)     = 0.0d0
      Amat(D1_RB)     = 0.0d0
      Amat2(D2_RA_RA) = 0.0d0
      Amat2(D2_RA_RB) = 0.0d0
      Amat2(D2_RB_RB) = 0.0d0
c
      call nwpwxc_c_pw91lda_d2(1.0d-20,2,1,1.0d0,rho,ec,Amat,Amat2)
c
      FCLDA(_FXC_E)    = ec
      FCLDA(_FXC_RA)   = Amat(D1_RA)
      FCLDA(_FXC_RB)   = Amat(D1_RB)
      FCLDA(_FXC_RARA) = Amat2(D2_RA_RA)
      FCLDA(_FXC_RARB) = Amat2(D2_RA_RB)
      FCLDA(_FXC_RBRB) = Amat2(D2_RB_RB)
c
      end
C>
C> @}
