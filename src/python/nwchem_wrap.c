#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#include <Python.h>
#include <abstract.h>

#include "rtdb.h"
#include "macdecls.h"
#include "global.h"

static PyObject *NwchemError;

#define MAX_VAR_NAME    256

static int rtdb_handle;            /* handle to the rtdb */

static PyObject *
wrap_rtdb_open(PyObject *self, PyObject *args)
{
   const char *filename, *mode;
   int  handle;

   if (PyArg_Parse(args, "(ss)", &filename, &mode)) {
       if (!rtdb_open(filename, mode, &rtdb_handle)) {
	   PyErr_SetString(NwchemError, "rtdb_open failed");
	   return NULL;
       }
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_open(filename, mode)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_rtdb_close(PyObject *self, PyObject *args)
{
   const char *mode;
   int  handle, result;

   if (PyArg_Parse(args, "s", &mode)) {
       if (!(result = rtdb_close(rtdb_handle, mode))) {
	   PyErr_SetString(NwchemError, "rtdb_close failed");
	   return NULL;
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: rtdb_close(mode)");
       return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_pass_handle(PyObject *self, PyObject *args)
{
   if (!(PyArg_Parse(args, "i", &rtdb_handle))) {
      PyErr_SetString(PyExc_TypeError, "Usage: pass_handle(rtdb_handle)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}

static PyObject *wrap_rtdb_print(PyObject *self, PyObject *args)
{
   int handle, flag;

   if (PyArg_Parse(args, "i", &flag)) {
      if (!rtdb_print(rtdb_handle, flag)) 
	   PyErr_SetString(NwchemError, "rtdb_print failed");
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: rtdb_print(flag)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}


static PyObject *wrap_rtdb_put(PyObject *self, PyObject *args)
{
    int handle, i, list, list_len;
    int ma_type = -1;
    char *name;
    int* int_array;
    double *dbl_array;
    char *char_array, *pchar;
    void *array = 0;
    PyObject *obj, *option_obj;

    if ((PyTuple_Size(args) == 2) || (PyTuple_Size(args) == 3)) {
	obj = PyTuple_GetItem(args, 0);      /* get var name */
	PyArg_Parse(obj, "s", &name);
	obj = PyTuple_GetItem(args, 1);      /* get an array or single value */
	
	if (PyList_Check(obj)) 
	    list = 1; 
	else 
	    list = 0;

	if (list) {
	    list_len = PyList_Size(obj);
	    if (   PyInt_Check(PyList_GetItem(obj, 0)))  ma_type = MT_F_INT;
	    if ( PyFloat_Check(PyList_GetItem(obj, 0)))  ma_type = MT_F_DBL;
	    if (PyString_Check(PyList_GetItem(obj, 0)))  ma_type = MT_CHAR;
	} else {
	    list_len = 1;
	    if (   PyInt_Check(obj))  ma_type = MT_F_INT;
	    if ( PyFloat_Check(obj))  ma_type = MT_F_DBL;
	    if (PyString_Check(obj))  ma_type = MT_CHAR; 
	} 
      
	if (PyTuple_Size(args) == 3) {
	    option_obj = PyTuple_GetItem(args, 2);      /* get optional type */
	    if (!(PyArg_Parse(option_obj, "i", &ma_type))) {
		PyErr_SetString(PyExc_TypeError, 
				"Usage: rtdb_put(value or values,[optional type])");
		return NULL;
	    }
	}
	
	if (ma_type != MT_CHAR) {
	    if (!(array = malloc(MA_sizeof(ma_type, list_len, MT_CHAR)))) {
		PyErr_SetString(PyExc_MemoryError,
				"rtdb_put failed allocating work array");
		return NULL;
	    }
	}
	    
	switch (ma_type) {
	case MT_INT:
	case MT_F_INT:  
	case MT_BASE + 11:	/* Logical */
	    int_array = array;
	    for (i = 0; i < list_len; i++) {
		if (list) 
		    PyArg_Parse(PyList_GetItem(obj, i), "i", int_array+i);
		else 
		    PyArg_Parse(obj, "i", int_array+i);
	    }
	    break;

	case MT_DBL:  
	case MT_F_DBL:
	    dbl_array = array;
	    for (i = 0; i < list_len; i++) {
		if (list) 
		    PyArg_Parse(PyList_GetItem(obj, i), "d", dbl_array+i);
		else 
		    PyArg_Parse(obj, "d", dbl_array+i);
	    }
	    break;

         case MT_CHAR: 
	     if (list) 
		 PyArg_Parse(PyList_GetItem(obj, 0), "s", &char_array); 
	     else 
		 PyArg_Parse(obj, "s", &char_array); 
	     if (!(array = strdup(char_array))) {
		 PyErr_SetString(PyExc_MemoryError,
				"rtdb_put failed copying string");
		 return NULL;
	     }		 
	     list_len = strlen(array) + 1;
	     break;
	     
	default:
	    PyErr_SetString(NwchemError, "rtdb_put: ma_type is incorrect");
	    if (array) free(array);
	    return NULL;
	    break;
	}                
                    
	if (!(rtdb_put(rtdb_handle, name, ma_type, list_len, array))) {
	    PyErr_SetString(NwchemError, "rtdb_put failed");
	    if (array) free(array);
	    return NULL;
	}

    } else {
	PyErr_SetString(PyExc_TypeError, 
			"Usage: rtdb_put(value or values,[optional type])");
	if (array) free(array);
	return NULL;
    }
    Py_INCREF(Py_None);
    if (array) free(array);
    return Py_None;
}

PyObject *wrap_rtdb_get(PyObject *self, PyObject *args)
{
   int i, result;
   int nelem, ma_type;
   int *int_array;
   double *dbl_array;
   char *char_array, *pchar;
   char *name;
   char *format_str=0, date[26], format_char;
   PyObject *returnObj;
   void *array=0;
   int ma_handle, ind;

   if (PyArg_Parse(args, "s", &name)) {
       if (!rtdb_ma_get(rtdb_handle, name, &ma_type, &nelem, &ma_handle)) {
	   PyErr_SetString(NwchemError, "rtdb_ma_get failed");
	   return NULL;
       }
       if (!MA_get_pointer(ma_handle, &array)) {
	   PyErr_SetString(NwchemError, "rtdb_ma_get failed");
	   return NULL;
       }
       /*printf("name=%s ma_type=%d nelem=%d ptr=%x\n",name, ma_type, 
	      nelem, array);*/

       switch (ma_type) {
       case MT_F_INT:
       case MT_INT  : 
       case MT_BASE + 11  : 
	   format_char = 'i'; break;
       case MT_F_DBL: 
       case MT_DBL  : 
	   format_char = 'd'; break;
	   break;
       case MT_CHAR : 
	   format_char = 's'; break;
       default:
	   PyErr_SetString(NwchemError, "rtdb_get: ma type incorrect");
	   (void) MA_free_heap(ma_handle);
	   return NULL;
	   break;
       }
       
       if (!(format_str = malloc(nelem+3))) {
	   PyErr_SetString(PyExc_MemoryError,
			   "rtdb_get failed allocating format string");
	   (void) MA_free_heap(ma_handle);
	   return NULL;
       }

       ind = 0;
       if (nelem > 1) format_str[ind++] = '[';
       for (i = 0; i < nelem; i++, ind++)
	   format_str[ind] = format_char;
       if (nelem > 1) format_str[ind++] = ']';
       format_str[ind] = 0;

       switch (ma_type) {
       case MT_F_INT:
       case MT_INT  : 
       case MT_F_DBL: 
       case MT_DBL  : 
       case MT_BASE + 11  : 
	   returnObj = Py_VaBuildValue(format_str, array); break; 
       case MT_CHAR : 
	   returnObj = Py_BuildValue("s#", array, nelem-1); break;
       }
   }
   else {
       PyErr_SetString(PyExc_TypeError, "Usage: value = rtdb_get(name)");
       if (format_str) free(format_str);
       return NULL;
   }
   (void) MA_free_heap(ma_handle);
   if (format_str) free(format_str);
   return returnObj;
}

static PyObject *wrap_task_energy(PyObject *self, PyObject *args)
{
    char *filename, *mode;
    char *theory;
    double energy;
    PyObject *returnObj;
    
    if (PyArg_Parse(args, "s", &theory)) {
	if (!rtdb_put(rtdb_handle, "task:theory", MT_CHAR, 
		      strlen(theory)+1, theory)) {
	    PyErr_SetString(NwchemError, "task_energy: putting theory failed");
	    return NULL;
	}
	if (!task_energy_(&rtdb_handle)) {
	    PyErr_SetString(NwchemError, "task_energy: failed");
	    return NULL;
	}
	if (!rtdb_get(rtdb_handle, "task:energy", MT_F_DBL, 1, &energy)) {
	    PyErr_SetString(NwchemError, "task_energy: getting energy failed");
	    return NULL;
	}
    }
    else {
	PyErr_SetString(PyExc_TypeError, "Usage: task_energy(theory)");
	return NULL;
    }
    
    return Py_VaBuildValue("d", &energy);
}

static PyObject *wrap_ga_nodeid(PyObject *self, PyObject *args)
{
    int nodeid = ga_nodeid_();
    if (args) {
	PyErr_SetString(PyExc_TypeError, "Usage: nodeid()");
	return NULL;
    }

    return Py_VaBuildValue("i", &nodeid);
}
    

static PyObject *wrap_nw_inp_from_string(PyObject *self, PyObject *args)
{
   char *pchar;

   if (PyArg_Parse(args, "s", &pchar)) {
       if (!nw_inp_from_string(rtdb_handle, pchar)) {
	   PyErr_SetString(NwchemError, "input_parse failed");
	   return NULL;
      }
   }
   else {
      PyErr_SetString(PyExc_TypeError, "Usage: input_parse(string)");
      return NULL;
   }
   Py_INCREF(Py_None);
   return Py_None;
}


/******************************************************************************/
/******************************************************************************/
/******************************************************************************/


static struct PyMethodDef nwchem_methods[] = {
   {"rtdb_open",       wrap_rtdb_open, 0}, 
   {"rtdb_close",      wrap_rtdb_close, 0}, 
   {"pass_handle",     wrap_pass_handle, 0}, 
   {"rtdb_print",      wrap_rtdb_print, 0}, 
   {"rtdb_put",        wrap_rtdb_put, 0}, 
   {"rtdb_get",        wrap_rtdb_get, 0}, 
   {"task_energy",     wrap_task_energy, 0}, 
   {"input_parse",     wrap_nw_inp_from_string, 0}, 
   {"ga_nodeid",       wrap_ga_nodeid, 0}, 
   {NULL, NULL}
};

void initnwchem()
{
    PyObject *m, *d;
    m = Py_InitModule("nwchem", nwchem_methods);
    d = PyModule_GetDict(m);
    NwchemError = PyErr_NewException("nwchem.error", NULL, NULL);
    PyDict_SetItemString(d, "error", NwchemError);
}

