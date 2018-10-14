#ifndef eval_nvmix_integral_h
#define eval_nvmix_integral_h

#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>


double eval_nvmix_integral_c(double *lower, double *upper, double *U, int n,
			     int d, double *cholScale, double ZERO, double ONE);
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP U, SEXP n, SEXP d,
			 SEXP cholScale, SEXP ZERO, SEXP ONE);

#endif /* eval_nvmix_integral_h */
