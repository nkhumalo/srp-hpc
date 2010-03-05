#if (defined(CRAY) || defined(ARDENT) || defined(WIN32))&& !defined(__crayx1)
#   define ga_access_callback_release_ GA_ACCESS_CALLBACK_RELEASE 
#endif

#include "ga.h"
#include "macdecls.h"
#if defined(CRAY)&& !defined(__crayx1) 
#define FATR
#include <fortran.h> /* Required for Fortran-C string interface on Crays */
#endif /* CRAY */
#ifdef WIN32
#include "typesf2c.h"
#endif
			     
/*\ PROVIDE ACCESS TO A PATCH OF A GLOBAL ARRAY WITH CALLBACK AND RELEASE
\*/
void FATR ga_access_callback_release_(g_a, ilo, ihi, jlo, jhi, 
				     callback, 
				     arg1, arg2, arg3, arg4, arg5, arg6, arg7)
     Integer *g_a, *ilo, *ihi, *jlo, *jhi;
     Integer (*callback)(Integer *,Integer *,Integer *,Integer *,Integer *,
			 void *, Integer*, 
			 void *, void *, void *, void *, void *, void *, void *);
     void *arg1, *arg2, *arg3, *arg4, *arg5, *arg6, *arg7;
{
  Integer lo[2], hi[2],ndim=GA_Ndim(*g_a), ld[2], result;
  int alo[2], ahi[2];
  void *ptr;

  if(ndim != 2) 
    GA_Error("ga_access: 2D API cannot be used for array dimension",ndim);

  lo[0]=*ilo;
  lo[1]=*jlo;
  hi[0]=*ihi;
  hi[1]=*jhi;
  nga_access_ptr(g_a,lo,hi,&ptr,ld);
  result = callback(g_a,ilo,ihi,jlo,jhi,ptr,&ld[0],arg1,arg2,arg3,arg4,arg5,arg6,arg7);
  alo[0]=lo[0];
  alo[1]=lo[1];
  ahi[0]=hi[0];
  ahi[1]=hi[1];
  if (result) {
    NGA_Release_update(*g_a, alo, ahi);
  } else {
    NGA_Release(*g_a, alo, ahi);
  }
} 

