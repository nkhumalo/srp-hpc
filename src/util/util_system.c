/*
 $Id: util_system.c,v 1.10 2005-01-06 02:10:40 edo Exp $
 */

#include <stdio.h>
#ifndef WIN32
#include <unistd.h>
#endif

extern int system(const char *);

#if defined(CRAY) && !defined(__crayx1)
#include <fortran.h>
#define FATR
#endif
#ifdef WIN32
#include "typesf2c.h"
#endif
#if (defined (LINUX) || defined (DECOSF)) && !defined(PSCALE)
#define __USE_BSD
#include <signal.h>
#undef __USE_BSD
#endif

typedef long Integer;		/*  FORTRAN integer */

#if defined(USE_FCD)
int fortchar_to_string(_fcd, int, char *, const int);
#else
int fortchar_to_string(const char *, int, char *, const int);
#endif

void ga_error(const char *, long);

#if defined(USE_FCD)
Integer FATR UTIL_SYSTEM(_fcd input)
{
    int lin  = _fcdlen(input);
#else
Integer util_system_(const char *input, int lin)
{
#endif
#if (defined(LINUX) || defined(DECOSF)) && !defined(PSCALE)
    int i;
	void (*Siginit)();
#endif
    char in[1024];
    if (!fortchar_to_string(input, lin, in, sizeof(in)))
	ga_error("util_system: fortchar_to_string failed for in",0);

#if (defined(LINUX) || defined(DECOSF)) && !defined(PSCALE)
	Siginit = signal(SIGCHLD,SIG_IGN);
	i = system(in);
	Siginit = signal(SIGCHLD,Siginit);
    return i;
#else
    return system(in);
#endif
}
