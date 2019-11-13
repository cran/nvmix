/* Register routines with R ***************************************************/

#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#include "eval_nvmix_integral.h"
#include "eval_densmix_integrand.h"
#include "logsumexp2.h"

static const R_CallMethodDef callMethods[] = {
    {"eval_nvmix_integral", (DL_FUNC) &eval_nvmix_integral, 10},
    {"eval_densmix_integrand", (DL_FUNC) &eval_densmix_integrand, 7},
    {"logsumexp2", (DL_FUNC) &logsumexp2, 3},
    {NULL, NULL, 0}
};

void R_init_nvmix(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, callMethods, NULL, NULL); /* s. WRE (2015, Section 5.4) */
    R_useDynamicSymbols(dll, FALSE);
}
