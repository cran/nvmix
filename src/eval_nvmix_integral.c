#include "eval_nvmix_integral.h"


/**
 * @title Evaluate Integrand for a Normal Variance Mixture
 * @param lower d-vector of lower evaluation limits
 * @param upper d-vector of upper evaluation limits
 * @param U vector representing '(n, d+1)-matrix' of uniforms
 *        (e.g. Sobol pointset) and realizations of sqrt(mix).
 *        The '(n, d+1)-matrix' is of the form [M1, U, M2], where the first
 *        column (M1) consists of realizations of sqrt(mix), the last column
 *        (M2) consists of the antithetic realizations of (M1) and the middle
 *        part is an (n, d-1)-matrix of uniforms (e.g. a Sobol point-set).
 * @param n sample size (i.e., number of rows of U)
 * @param d dimension minus 1 (i.e., number of columns of U minus 1); this is
 *        the dimension of the original problem
 * @param factor lower triangular Cholesky factor as vector
 * @param ZERO smallest number x > 0 such that x != 0
 * @param ONE   largest number x < 1 such that x != 1
 * @return mean estimate mean(f(U)) of E(f(U)) using antithetic variates
 * @note lower and upper can have -/+Inf entries. While pnorm() would give the
 correct result, it safes time to check each time if the argument is -/+Inf
 *       and setting the value to 0/1 rather than calling pnorm().
 * @author Erik Hintz and Marius Hofert
 */
double eval_nvmix_integral_c(double *lower, double *upper, double *U, int n, int d,
                             double *factor, double ZERO, double ONE)
{
    double yorg[d-1], sqrtmixorg, dorg, difforg, forg, scprodorg;
    double yant[d-1], sqrtmixant, dant, diffant, fant, scprodant;
    /* Note: <name>org stands for "original", <name>ant for antithetic */
    /* y:       vector to save phi^{-1}(dj+uj(ej-dj)) */
    /* sqrtmix: used to store sqrt(F_w^{-1}(u_0)) */
    /* d:       current values of di from the paper */
    /* diff:    current values of (ei-di) from the paper */
    /* f:       current value of (e1-d1) * (e2-d2) * ... * (ei-di) */
    /* scprod:  scalar product sum factor_{ij} y_j */

    double tmp; /* to store temporary values */
    double mean = 0; /* to store the result */
    int i, j, l; /* counters for loops */

    /* For each row of U (for access, use U[s,k] = U[k * numrows + s]) */
    for(j = 0; j < n; j++){

        /* Grab the realizations of sqrt(mix) = sqrt(W) */
        sqrtmixorg = U[j];
        sqrtmixant = U[d*n+j];

        /* Check if entry of lower is -Inf */
        if(lower[0] == R_NegInf){
            dorg = 0;
            dant = 0;
        } else {
            dorg = pnorm(lower[0] / (factor[0] * sqrtmixorg), 0, 1, 1, 0);
            dant = pnorm(lower[0] / (factor[0] * sqrtmixant), 0, 1, 1, 0);
        }

        /* Check if entry of b is +Inf */
        if(upper[0] == R_PosInf){
            difforg = 1 - dorg;
            diffant = 1 - dant;
        } else {
            difforg = pnorm(upper[0] / (factor[0] * sqrtmixorg), 0, 1, 1, 0) - dorg;
            diffant = pnorm(upper[0] / (factor[0] * sqrtmixant), 0, 1, 1, 0) - dant;
        }

        /* Go through all d-1 columns (without first and last) */
        /* For better readability, we start at i = 0 */
        forg = difforg;
        fant = diffant;
        for(i = 0; i < d-1; i++){
            /* U[i * n + j] corresponds to U[j,i] in the orginal matrix */
            tmp = dorg + U[(i+1) * n + j] * difforg;

            /* Check if too close to 0 or 1 */
            if(tmp < ZERO){
                tmp = ZERO;
            }
            if(tmp > ONE){
                tmp = ONE;
            }
            yorg[i] = qnorm(tmp, 0, 1, 1, 0);

            /* The same for the antithetic value */
            tmp = dant + (1-U[(i+1) * n + j]) * diffant;
            if(tmp < ZERO){
                tmp = ZERO;
            }
            if(tmp > ONE){
                tmp = ONE;
            }
            yant[i] = qnorm(tmp, 0, 1, 1, 0);

            /* Calculate the scalar product sum factor[i,j] y[j] for j = 1:(i-1) */
            scprodorg = 0;
            scprodant = 0;
            for(l = 0; l < (i+1); l++){
                /* factor[l * d + i+1] corresponds to factor[i+1,l] in the original matrix */
                scprodorg += yorg[l] * factor[l * d + i+1];
                scprodant += yant[l] * factor[l * d + i+1];
            }

            /* Calculate new d and diff = e-d */
            if(lower[i+1] == R_NegInf){
                dorg = 0;
                dant = 0;
            } else {
                dorg = pnorm((lower[i+1] / sqrtmixorg - scprodorg) / factor[(i+1)*(d+1)], 0, 1, 1, 0);
                dant = pnorm((lower[i+1] / sqrtmixant - scprodant) / factor[(i+1)*(d+1)], 0, 1, 1, 0);
            }
            if(upper[i+1] == R_PosInf){
                difforg = 1 - dorg;
                diffant = 1 - dant;
            } else {
                difforg = pnorm((upper[i+1] / sqrtmixorg - scprodorg) / factor[(i+1)*(d+1)], 0, 1, 1, 0) - dorg;
                diffant = pnorm((upper[i+1] / sqrtmixant - scprodant) / factor[(i+1)*(d+1)], 0, 1, 1, 0) - dant;
            }
            forg *= difforg;
            fant *= diffant;
        }
        mean += (forg+fant)/2;
    }
    mean = mean / n;
    return(mean);
}

/**
 * @title R Interface for eval_nvmix_integral_c()
 * @param see eval_nvmix_integral_c()
 * @return mean(f(U)) where f is the integrand and U specifies the point-set
 * @author Erik Hintz, Marius Hofert (polishing)
 */
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP U, SEXP n, SEXP d,
			 SEXP factor, SEXP ZERO, SEXP ONE)
{
    double res = eval_nvmix_integral_c(REAL(lower), REAL(upper), REAL(U),
				       INTEGER(n)[0], INTEGER(d)[0],
				       REAL(factor), REAL(ZERO)[0],
                                       REAL(ONE)[0]);
    return ScalarReal(res);
}
