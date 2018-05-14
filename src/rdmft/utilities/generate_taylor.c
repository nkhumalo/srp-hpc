/*----------------------------------------------------------------------

Author: Hubertus van Dam    Date: 04/27/2018

Generate Taylor coefficients of ln(x!)+ln((1-x)!)

This program uses the GNU Scientific Library which is licensed under
GNU GPL. Therefore this program is also licenced under the GNU
Public License.

This program generates Taylor coefficients of the function

  f(x) = ln(x!)+ln((1-x)!)

where 0 <= x <= 1. The Taylor series in question is expanded around
a=1/4. The program generates coefficients 0 through to 28.
----------------------------------------------------------------------*/
#include "stdio.h"
#include "gsl/gsl_sf_gamma.h"

void main(void)
{
    int    mxcoeff = 28;
    double coeff;
    double error;
    double a   = 0.25e0;
    double one = 1.00e0;
    double two = 2.00e0;
    double facinvb = 1.0e0; // bottom range i.e. h=-0.25
    double facinvt = 1.0e0; // top    range i.e. h=+0.25
    int    ii;           // counter
    int    ierr;
    gsl_sf_result out;
    printf("         Taylor coeff            ! term  coeff error     est term size  est term size\n");
    for (ii = 0; ii <= mxcoeff; ii++) {
        if (ii==0) {
            coeff = gsl_sf_lngamma(one+a)
                  + gsl_sf_lngamma(two-a);
        }
        else {
            facinvb /= (1.0e0*ii);
            facinvb *= -0.25e0;
            facinvt /= (1.0e0*ii);
            facinvt *= +0.25e0;
            if (ierr = gsl_sf_psi_n_e(ii-1,one+a,&out)) {
                printf("error a: %d %24.16e %24.16e \n",
                       ierr,out.val,out.err);
            }
            else {
                coeff = out.val;
                error = out.err;
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
                error += out.err;
            }
        }
        printf("     & %24.16e, ! %3d %14.6e   %14.6e %14.6e\n",coeff,ii,error,coeff*facinvb,coeff*facinvt);
    }
}

//----------------------------------------------------------------------
