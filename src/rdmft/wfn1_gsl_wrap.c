#ifdef USE_GSL
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#endif
#include <stdlib.h>
#include <stdio.h>

double gsl_sf_lngamma_(double* x)
{
#ifdef USE_GSL
   return gsl_sf_lngamma(*x);
#else
   printf("gsl_sf_lngamma_: USE_GSL not set\n");
   exit(1);
#endif
}

double gsl_sf_psi_(double* x)
{
#ifdef USE_GSL
   return gsl_sf_psi(*x);
#else
   printf("gsl_sf_psi_: USE_GSL not set\n");
   exit(1);
#endif
}

double gsl_sf_gamma_(double* x)
{
#ifdef USE_GSL
   return gsl_sf_gamma(*x);
#else
   printf("gsl_sf_gamma_: USE_GSL not set\n");
   exit(1);
#endif
}
