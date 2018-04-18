#ifdef USE_GSL
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>
#endif

double gsl_sf_lngamma_(double* x)
{
#ifdef USE_GSL
   return gsl_sf_lngamma(*x);
#else
   exit(1);
#endif
}

double gsl_sf_psi_(double* x)
{
#ifdef USE_GSL
   return gsl_sf_psi(*x);
#else
   exit(1);
#endif
}
