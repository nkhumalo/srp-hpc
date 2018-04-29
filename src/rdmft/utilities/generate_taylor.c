/*----------------------------------------------------------------------

Author: Hubertus van Dam    Date: 04/27/2018

Generate Taylor coefficients of ln(x!)+ln((1-x)!)

This program uses the GNU Scientific Library which is licensed under
GNU GPL. Therefore this program is also licenced under the GNU
Public License.

This program generates Taylor coefficients of the function

  f(x) = ln(x!)+ln((1-x)!)

where 0 <= x <= 1. The Taylor series in question is expanded around
a=1/4. The program generates coefficients 0 through to 27.
----------------------------------------------------------------------*/
#include "stdio.h"
#include "gsl/gsl_sf_gamma.h"

void main(void)
{
    int    mxcoeff = 27;
    double coeff;
    double a   = 0.25e0;
    double one = 1.00e0;
    double two = 2.00e0;
    int    ii;           // counter
    int    ierr;
    gsl_sf_result out;
    for (ii = 0; ii <= mxcoeff; ii++) {
        if (ii==0) {
            coeff = gsl_sf_lngamma(one+a)
                  + gsl_sf_lngamma(two-a);
        }
        else {
            if (ierr = gsl_sf_psi_n_e(ii-1,one+a,&out)) {
                printf("error a: %d %24.16e %24.16e \n",
                       ierr,out.val,out.err);
            }
            else {
                coeff = out.val;
            }
            if (ierr = gsl_sf_psi_n_e(ii-1,two-a,&out)) {
                printf("error b: %d %24.16e %24.16e \n",
                       ierr,out.val,out.err);
            }
            else {
                if (ii % 2 == 0) {
                    coeff += out.val;
                }
                else {
                    coeff -= out.val;
                }
            }
        }
        printf("     & %24.16e, ! %3d\n",coeff,ii);
    }
}

//----------------------------------------------------------------------
