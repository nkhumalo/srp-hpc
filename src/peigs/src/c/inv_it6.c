/*
 $Id: inv_it6.c,v 1.3 1999-07-28 00:39:25 d3e129 Exp $
 *======================================================================
 *
 * DISCLAIMER
 *
 * This material was prepared as an account of work sponsored by an
 * agency of the United States Government.  Neither the United States
 * Government nor the United States Department of Energy, nor Battelle,
 * nor any of their employees, MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR
 * ASSUMES ANY LEGAL LIABILITY OR RESPONSIBILITY FOR THE ACCURACY,
 * COMPLETENESS, OR USEFULNESS OF ANY INFORMATION, APPARATUS, PRODUCT,
 * SOFTWARE, OR PROCESS DISCLOSED, OR REPRESENTS THAT ITS USE WOULD NOT
 * INFRINGE PRIVATELY OWNED RIGHTS.
 *
 * ACKNOWLEDGMENT
 *
 * This software and its documentation were produced with Government
 * support under Contract Number DE-AC06-76RLO-1830 awarded by the United
 * States Department of Energy.  The Government retains a paid-up
 * non-exclusive, irrevocable worldwide license to reproduce, prepare
 * derivative works, perform publicly and display publicly by or for the
 * Government, including the right to distribute to other Government
 * contractors.
 *
 *======================================================================
 *
 *  -- PEIGS  routine (version 2.1) --
 *     Pacific Northwest Laboratory
 *     July 28, 1995
 *
 *======================================================================
 */
/* **********************************************
   
   PeIGS internal routine
   
   Integer inv_it4( n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, eval, eps, stpcrt, onenrm, iwork, work)
   Integer *n, *c1, *cn, *b1, *bn, *Zbegin, *map, *mapvec, *iwork;
   DoublePrecision *d, *e, **vector, *eval, *work, *eps, *stpcrt, *onenrm;
   
   This routine performs inverse iteration.
   
   */   

#include <stdio.h>
#include <math.h>

#include "globalp.c.h"

#include "clustr_inv.h"

#define    ITS   1
#define    MAXITR 100

Integer inv_it6( n, c1, cn, b1, bn, Zbegin, map, mapvec, vector, d, e, eval, eps, stpcrt, onenrm, iwork, work)
     Integer *n, *c1, *cn, *b1, *bn, *Zbegin, *map, *mapvec, *iwork;
     DoublePrecision *d, *e, **vector, *eval, *work, *eps, *stpcrt, *onenrm;
     
     /*
       this routine performs inverse iteration on *n vectors given by map[0:*n-1]
       and whose storage location on each processor is given by mm
       
       n = dimension of the vectors
       m = number of vectors to inverse iterate
       
       d is the diag of the matrix
       e is the sub and super diagonals of the matrix
       eval is the list of the eigenvalues  eval(c1) ... eval[cn] are the eigenvalues; mapvec[k] gives
       the true index of the eigenvalues
       
       */
{
  static Integer IONE = 1, IMINUSONE = -1;
  Integer j, blksz, i_1, niter;
  Integer nrmchk, indrv1, indrv2, indrv3, indrv4, indrv5;
  Integer jmax, k, lsiz, me;
  DoublePrecision xj, xjm, sep, nrm, scl;
  DoublePrecision tol, *ptr;
  Integer info, mxmynd();
  Integer csiz, ibad;
  Integer indx22, i, lll;
  DoublePrecision delta = 0., dummy;
  DoublePrecision ztz, *workk, *ld, *lld, *dptr, *wwork1;
  Integer zbegin1, kk, zend, bb1 = *b1+1, bbn = *bn+1, *iwork1;
  
  extern void dscal_(), dlagtf_(), dlagts_(), dcopy_(), dscal_();
  extern Integer idamax_();
  extern DoublePrecision dnrm2_(), dasum_(), ddot_();
  extern Integer mxmynd_();
  extern DoublePrecision dgetavec_();
  
  
  extern void mgs_prev ();
  
  me = mxmynd_();
  
  blksz = *bn - *b1 + 1;
  
  lsiz = *n;
  csiz = *cn - *c1 + 1;
  
  iwork1 = (Integer *) malloc( 16* lsiz *sizeof(Integer));
  wwork1 = (DoublePrecision *) malloc( 16* lsiz *sizeof(DoublePrecision));
  
  for(i = 0;i < 16*lsiz;i++){
    iwork1[i] = 0;
    wwork1[i] = 0.;
  }
  
  indrv1 = 0;
  indrv2 = lsiz + indrv1;
  indrv3 = lsiz + indrv2;
  indrv4 = lsiz + indrv3;
  indrv5 = lsiz + indrv4;
  
  ld = wwork1 + 10* lsiz;
  lld = wwork1 + 11* lsiz;
  
  ibad = 0;

#ifdef DEBUG  
  fprintf(stderr, " csiz = %d blksize = %d \n", csiz, blksz);
#endif

  for(i = *b1;i < *bn+1;i++){
    ld[i] = d[i]*e[i];
    lld[i] = ld[i]*e[i];
  }
  /*
    iwork1 = iwork+*n;
    */
  
  sep = 0.;
  k = *Zbegin;  
  i_1 = *c1;
  for ( j = 0; j < csiz ;  j++ ) {
    if ( map[ i_1 ] == me ) {
      mapvec[k] = i_1;
      xj = eval[i_1];
      dgetavec2_( &j, &xj, &delta, n, &bb1, &bbn, e, d, ld, lld, wwork1,
		  wwork1 + lsiz, wwork1 + 2*lsiz, wwork1 + 3*lsiz, wwork1 + 4*lsiz,
		  wwork1 + 5* lsiz, wwork1 + 6* lsiz,
		  vector[k] + *b1, &kk, &ztz, &zbegin1, &zend, &j, iwork1 );
      printf(" invit6 me = %d csiz = %d j = %d xj = %f c1 = %d cn = %d ztz = %f \n", me,
	     csiz, j, xj, *c1, *cn, ztz);
      ztz = 1./ztz;
      
      dscal_(&blksz, &ztz, vector[k] + *b1, &IONE);
      i_1++;
      k++;
    }
  }
  
  free(wwork1);
  free(iwork1);
  return(ibad);
}  
