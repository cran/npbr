#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void sigma2m(void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"sigma2m", (DL_FUNC) &sigma2m, 6},
    {NULL, NULL, 0}
};

void R_init_npbr(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

