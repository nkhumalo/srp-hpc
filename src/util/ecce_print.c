#include <stdio.h>
#include <string.h>
#include "types.f2c.h"
#include "ecce_print.h"
#include "macommon.h"

Integer MA_sizeof(Integer, Integer, Integer);

/* These two routines defined in rtdb_seq.c */
void ma_print(FILE *file, const int ma_type, const int nelem, const void *p);
const char *ma_typename(const int ma_type);

static int ecce_print_enabled;	/* If true then print */

static FILE *ecce_file;		/* ECCE output file  */

static char module_stack[4096];

static void module_stack_push(const char *module)
/*
  Push the module onto the stack.  Module is assumed to not contain
  white space.  

  If the push fails print a warning message once only.
  */
{
    int stacklen = strlen(module_stack);
    int modlen = strlen(module);
    static int print_warning = 1;
    
    if ( (modlen+stacklen+2) >= sizeof(module_stack) ) {
	if (print_warning) {
	    fprintf(stderr, "!! ecce_print: module_stack too small %d %d\n",
		    modlen+stacklen+2, sizeof(module_stack));
	    print_warning = 0;
	}
    }
    else {
	if (stacklen) {
	    (void) strcpy(module_stack+stacklen, " ");
	    stacklen++;
	}
	(void) strcpy(module_stack+stacklen, module);
    }
}

static void module_stack_pop(const char *module)
/*
  Pop the stack at module name.  Print at most one message per
  calculation if module is not on the stack.
  */
{
    int stacklen;
    int modlen = strlen(module);
    static int print_warning = 1;

    while ( (stacklen = strlen(module_stack)) ) {
	char *m = module_stack+stacklen-1;
	int foundit;

	while ( (m > module_stack) && (*m != ' '))
	    m--;		/* Find beginning of last name on stack */
	if (*m == ' ') m++;
	
	foundit = (strncmp(m, module, modlen) == 0);

	if (m > module_stack) m--; /* To remove preceeding space */
	*m = 0;

	if (foundit) return;
    }

    if (print_warning) {
	fflush(stdout);
	fprintf(stderr," !! ecce_print: module (%s) is not on stack (%s)\n",
		module, module_stack);
	print_warning = 0;
    }

    module_stack[0] = 0;	/* Not found ... pop entire stack since this
				   will eventually yield the correct state */
}

static void remove_blanks(char *t)
{
    while (*t) {
	if (*t == ' ') *t = '_';
	t++;
    }
}

static void print_info(const char *info, const char *key, const char *type,
		       int dim1, int dim2)
/*
  Print

  <module stack>%<info>%<key>%dim1 ... %<type>

  If (dim2 <= 0) do not print it.  If more than two dims are needed
  we should probably use varargs but that's too much like hard work for now.
*/
{
    fprintf(ecce_file, "%s%%%s%%%s%%%d", module_stack, info, key, dim1);
    if (dim2) fprintf(ecce_file, " %d", dim2);
    fprintf(ecce_file, "%%%s\n",type);
}
		       

void ecce_print_module_entry(const char *module) 
{
    char buf[256];

    if (!ecce_print_enabled) return;

    strncpy(buf, module, 255); buf[255] = 0; /* Quietely truncate long names */
    remove_blanks(buf);
    
    print_info("begin", "entry", "char", 1, 0);
    fprintf(ecce_file, "%s\n", buf);
    print_info("end", "entry", "char", 1, 0);
    fflush(ecce_file);

    module_stack_push(buf);
}

void ecce_print_module_exit(const char *module, const char *status) 
{
    char buf[256];

    if (!ecce_print_enabled) return;

    strncpy(buf, module, 255); buf[255] = 0; /* Quietly truncate long names */
    remove_blanks(buf);
    
    print_info("begin", "exit", "char", 2, 0);
    fprintf(ecce_file, "%s\n%s\n", buf, status);
    print_info("end", "exit", "char", 1, 0);
    fflush(ecce_file);

    module_stack_pop(buf);
}
    
void ecce_print1(const char *key, int ma_type, const void *data, int dim1)
{
    static int print_warning = 1;
    const char *typename = ma_typename(ma_type);
    char **c;
    int i;

    if (!ecce_print_enabled) return;
    if (strcmp(typename,"invalid") == 0) {
	if (print_warning) {
	    fprintf(stderr,"!! ecce_print: invalid type %d for %s\n", 
		    ma_type, key);
	    print_warning = 0;
	}
	return;
    }

    print_info("begin", key, typename, dim1, 0);
    if (strcmp(typename, "char"))
	ma_print(ecce_file, ma_type, dim1, data);
    else
	for (c=(char **) data, i=0; i<dim1; i++)
	    fprintf(ecce_file, "%s\n", c[i]);
    print_info("end", key, typename, dim1, 0);
    fflush(ecce_file);
}


void ecce_print2(const char *key, int ma_type, 
		 const void *data, int ld1, int dim1, int dim2)
{
    static int print_warning = 1;
    const char *typename = ma_typename(ma_type);
    const char *cdata = data;
    int typesize = MA_sizeof(ma_type, 1, MT_C_CHAR);
    int i;

    if (!ecce_print_enabled) return;
    if (strcmp(typename,"invalid") == 0) {
	if (print_warning) {
	    fprintf(stderr,"!! ecce_print: invalid type %d for %s\n", 
		    ma_type, key);
	    print_warning = 0;
	}
	return;
    }

    print_info("begin", key, typename, dim1, dim2);
    for (i=0; i<dim2; i++, cdata+= (ld1*typesize)) 
	ma_print(ecce_file, ma_type, dim1, (void *) cdata);
    print_info("end", key, typename, dim1, dim2);
    fflush(ecce_file);
}
    
void ecce_print_control(int new, int *old)
{
    *old = ecce_print_enabled;
    ecce_print_enabled = (new && ecce_file);
}
    
void ecce_print_file_open(const char *filename)
{
    if (!(ecce_file = fopen(filename, "w+"))) {
	fprintf(stderr,"!! ecce_print: failed to open %s\n", filename);
	ecce_print_enabled = 0;
	return;
    }
    ecce_print_enabled = 1;
}

void ecce_print_file_close(void)
{
    if (ecce_file) (void) fclose(ecce_file);
    ecce_file = (FILE *) 0;
    ecce_print_enabled = 0;
}


/*****************************************************************
  Following stuff is FORTRAN wrapping for the above
  ****************************************************************/


static int fortchar_to_string(const char *f, int flen, char *buf, 
			      const int buflen)
{
  while (flen-- && f[flen] == ' ')
    ;

  if ((flen+1) >= buflen)
    return 0;			/* Won't fit */

  flen++;
  buf[flen] = 0;
  while(flen--)
    buf[flen] = f[flen];

  return 1;
}

#if defined(CRAY) || defined(CRAY_T3D)
#define ecce_print_file_open_    ECCE_PRINT_FILE_OPEN
#define ecce_print_file_close_   ECCE_PRINT_FILE_CLOSE
#define ecce_print_control_      ECCE_PRINT_CONTROL
#define ecce_print1_             ECCE_PRINT1
#define ecce_print2_             ECCE_PRINT2
#define ecce_print1_char_        ECCE_PRINT1_CHAR
#define ecce_print_module_entry_ ECCE_PRINT_MODULE_ENTRY
#define ecce_print_module_exit_  ECCE_PRINT_MODULE_EXIT
#endif

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print_file_open_(_fcd f) 
{
    const char *filename = _fcdtocp(f);
    int flen = _fcdlen(f);
#else
void ecce_print_file_open_(const char *filename, int flen)
{
#endif
    char buf[1024];

    if (!fortchar_to_string(filename, flen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print_file_open: name too long? (%d %d)\n",
		flen, sizeof(buf));
	return;
    }

    ecce_print_file_open(buf);
}

void ecce_print_file_close_(void)
{
    ecce_print_file_close();
}

void ecce_print_control_(Integer *pnew, Integer *pold)
{
    int old;
    ecce_print_control((int) *pnew, &old);
    *pold = (Integer) old;
}

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print2_(fcd f, Integer *ma_type, 
		  const void *data, Integer *ld1, Integer *dim1, Integer *dim2)
{
    const char *key = _fcdtocp(f);
    int keylen = _fcdlen(f);
#else
void ecce_print2_(const char *key, Integer *ma_type, 
		  const void *data, Integer *ld1, Integer *dim1, Integer *dim2,
		  int keylen)
{
#endif
    char buf[1024];

    if (!fortchar_to_string(key, keylen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print: key too long (%d %d)\n", 
		keylen, sizeof(buf));
	return;
    }

    ecce_print2(buf, (int) *ma_type, data, (int) *ld1, (int) *dim1, 
		(int) *dim2);
}

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print1_(fcd f, Integer *ma_type, const void *data, Integer *dim1)
{
    const char *key = _fcdtocp(f);
    int keylen = _fcdlen(f);
#else
void ecce_print1_(const char *key, Integer *ma_type, 
		  const void *data, Integer *dim1, int keylen)
{
#endif
    char buf[1024];

    if (!fortchar_to_string(key, keylen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print: key too long (%d %d)\n", 
		keylen, sizeof(buf));
	return;
    }

    ecce_print1(buf, (int) *ma_type, data, (int) *dim1);
}

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print_module_entry_(_fcd f) 
{
    const char *module = _fcdtocp(f);
    int modlen = _fcdlen(f);
#else
void ecce_print_module_entry_(const char *module, int modlen) 
{
#endif
    char buf[1024];

    if (!fortchar_to_string(module, modlen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print_module_entry: name too long? (%d %d)\n",
		modlen, sizeof(buf));
	return;
    }

    ecce_print_module_entry(buf);
}

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print_module_exit_(_fcd f, _fcd g) 
{
    const char *module = _fcdtocp(f);
    int modlen = _fcdlen(f);
    const char *status = _fcdtocp(g);
    int statlen = _fcdlen(g);
#else
void ecce_print_module_exit_(const char *module, const char *status,
			     int modlen, int statlen) 
{
#endif
    char buf[1024], buf1[1024];

    if (!fortchar_to_string(module, modlen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print_module_exit: name too long? (%d %d)\n",
		modlen, sizeof(buf));
	return;
    }

    if (!fortchar_to_string(status, statlen, buf1, sizeof(buf1))) {
	fprintf(stderr,"!! ecce_print_module_exit: status too long? (%d %d)\n",
		statlen, sizeof(buf1));
	return;
    }

    ecce_print_module_exit(buf, buf1);
}

#if defined(CRAY) || defined(CRAY_T3D)
void ecce_print1_char_(fcd f, fcd g, Integer *dim1)
{
    const char *key = _fcdtocp(f);
    int keylen = _fcdlen(f);
    const char *data = _fdctocp(g);
    int dlen = _fcdlen(g);
#else
void ecce_print1_char_(const char *key, const char *data, Integer *dim1, 
		       int keylen, int dlen)
{
#endif
    char buf[1024], buf1[1024];
    int i;

    if (!fortchar_to_string(key, keylen, buf, sizeof(buf))) {
	fprintf(stderr,"!! ecce_print: key too long (%d %d)\n", 
		keylen, sizeof(buf));
	return;
    }

    print_info("begin", key, "char", (int) *dim1, 0);
    for (i=0; i<*dim1; i++, data+= dlen) {
	if (!fortchar_to_string(data, dlen, buf1, sizeof(buf1))) {
	    fprintf(stderr,"!! ecce_print1_char: datum too long (%d %d)\n", 
		    dlen, sizeof(buf1));
	    return;
	}
	fprintf(ecce_file,"%s\n", buf1);
    }
    print_info("end", key, "char", (int) *dim1, 0);
    fflush(ecce_file);
}

