#ifndef eval_nvmix_integral_h
#define eval_nvmix_integral_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


void eval_nvmix_integral_c(double *lower, double *upper, double *U, int n, int d,
                           int r, int *kfactor, double *factor, double ZERO, double ONE, double *res);
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP U, SEXP n, SEXP d, SEXP r, SEXP kfactor,
			 SEXP factor, SEXP ZERO, SEXP ONE);

#endif /* eval_nvmix_integral_h */


