#include <stdio.h>
#include "mdglobal.h"
#include "md5.h"
#if defined(CRAY)&& !defined(__crayx1)
#include <fortran.h>
#endif
#include "typesf2c.h"

/*
 $Id: md5wrap.c,v 1.7 2003-08-13 18:06:11 edo Exp $
 */

#if defined(USE_FCD)
extern int string_to_fortchar(_fcd f, int flen, char *buf);
#else
extern int string_to_fortchar(char *f, int flen, char *buf);
#endif

#if (defined(CRAY) || defined(WIN32))&& !defined(__crayx1)
#define checksum_simple_ CHECKSUM_SIMPLE
#define checksum_init_   CHECKSUM_INIT
#define checksum_update_ CHECKSUM_UPDATE
#define checksum_final_  CHECKSUM_FINAL
#define checksum_char_update_ CHECKSUM_CHAR_UPDATE
#define checksum_char_simple_ CHECKSUM_CHAR_SIMPLE
#endif

/* This passes mddriver -x test using munged MDstring in mddriver.c 
   and also produces the same checksum when compressing stdin */

/* Define this to compile the test main program */
/*#define MD5TEST*/

static MD5_CTX context;
#ifdef EXT_INT
typedef long integer;		/* FORTRAN integer */
#else
typedef int integer;		/* FORTRAN integer */
#endif

/* C interface */

void checksum_init(void)
{
    MD5Init(&context);
}

void checksum_update(int len, const void *buf)
{
    MD5Update(&context, buf, (unsigned int) len);
}

static void checksum_sum_to_string(const unsigned char sum[16], char csum[33])
{
    int i;
    char *c;

    for (c=csum, i=0; i < 16; i++, c+=2)
	(void) sprintf (c, "%02x", sum[i]);
    *c = 0;
}

void checksum_final(char csum[33])
{
    unsigned char sum[16];

    MD5Final(sum, &context);
    checksum_sum_to_string(sum, csum);
}

void checksum_simple(int len, const void *buf, char csum[33])
{
    checksum_init();
    checksum_update(len, buf);
    checksum_final(csum);
}

/* Don't need this routine ? */
static void checksum_string_to_sum(const char csum[33], unsigned char sum[16])
{
    int i, j;
    char buf[3];
    const char *c;

    buf[2] = 0;

    for (c=csum, i=0; i < 16; i++, c+=2) {
	buf[0] = c[0]; buf[1] = c[1];
	if (sscanf(buf, "%x", &j) != 1) {
	    fprintf(stderr,"checksum_sum_to_string: sscanf failed!\n");
	    exit(1);
	}
	sum[i] = j;
    }
}

/* Fortran interface */

void FATR checksum_init_(void)
{
    checksum_init();
}

void FATR checksum_update_(integer *len, const void *buf)
{
    checksum_update((int) *len, buf);
}

#if defined(USE_FCD)
void FATR checksum_char_update_(_fcd f)
{
    checksum_update(_fcdlen(f), _fcdtocp(f));
}
#else
void FATR checksum_char_update_(const char *buf, int len)
{
    checksum_update(len, buf);
}
#endif

#if defined(USE_FCD)
void FATR checksum_final_(_fcd f)
{
    int flen = _fcdlen(f);
#else
void FATR checksum_final_(char *f, int flen)
{
#endif
    char tmp[33];
    
    checksum_final(tmp);

    if (!string_to_fortchar(f, flen, tmp)) {
	fprintf(stderr,"checksum_final_: sum needs 32 chars have %d\n", flen);
	exit(1);
    }
}

#if defined(USE_FCD)
void checksum_simple_(integer *len, const void *buf, _fcd f)
{
    checksum_init();
    checksum_update((int) *len, buf);
    checksum_final_(f);
}
#else
void checksum_simple_(integer *len, const void *buf, char *f, int flen)
{
    checksum_init();
    checksum_update((int) *len, buf);
    checksum_final_(f, flen);
}
#endif

#if defined(USE_FCD)
void checksum_char_simple_(_fcd b, _fcd f)
{
    checksum_init();
    checksum_char_update_(b);
    checksum_final_(f);
}
#else
void checksum_char_simple_(const char *buf, char *sum, int blen, int slen)
{
    checksum_init();
    checksum_char_update_(buf, blen);
    checksum_final_(sum, slen);
}
#endif

#ifdef MD5TEST    
int cmain()
{
    char buf[1024], csum[33];
    int len;

    checksum_init();
    while ((len = fread(buf, 1, sizeof(buf), stdin)) > 0)
	checksum_update(len, buf);
    checksum_final(csum);

    printf("sum = %s\n", csum);

    return 0;
}

/* Normally defined by NWChem */

#if defined(USE_FCD)
static int string_to_fortchar(_fcd f, int flen, char *buf)
#else
static int string_to_fortchar(char *f, int flen, char *buf)
#endif
{
  int len = (int) strlen(buf), i;
#if defined(USE_FCD)
  flen = _fcdlen(f);
#endif

  if (len > flen) 
    return 0;			/* Won't fit */

#if defined(USE_FCD)
  for (i=0; i<len; i++)
    _fcdtocp(f)[i] = buf[i];
  for (i=len; i<flen; i++)
    _fcdtocp(f)[i] = ' ';
#else
  for (i=0; i<len; i++)
    f[i] = buf[i];
  for (i=len; i<flen; i++)
    f[i] = ' ';
#endif
  return 1;
}
#endif
