#include "eval_nvmix_integral.h"


/**
 * @title Evaluate Integrand for a Normal Variance Mixture
 * @param lower d-vector of lower evaluation limits
 * @param upper d-vector of upper evaluation limits
 * @param U vector representing '(n, r+1)-matrix' of uniforms
 *        (e.g. Sobol pointset) and realizations of sqrt(mix).
 *        The '(n, r+1)-matrix' is of the form [M1, U, M2], where the first
 *        column (M1) consists of realizations of sqrt(mix), the last column
 *        (M2) consists of the antithetic realizations of (M1) and the middle
 *        part is an (n, r-1)-matrix of uniforms (e.g. a Sobol point-set).
 * @param n sample size (i.e., number of rows of U)
 * @param d dimension of the original problem
 * @param r rank of factor
 * @param kfactor vector of length r giving height of each step in 'factor'
 * @param factor lower triangular Cholesky factor as vector
 * @param ZERO smallest number x > 0 such that x != 0
 * @param ONE   largest number x < 1 such that x != 1
 * @param res 2-vector to store results (mean + estimated variance)
 * @return 2-vector consisting of mean(f(U)) and var(f(U)) using antithetic variates
 * @note lower and upper can have -/+Inf entries. While pnorm() would give the
 correct result, it safes time to check each time if the argument is -/+Inf
 *       and setting the value to 0/1 rather than calling pnorm().
 * @author Erik Hintz and Marius Hofert
 */
void eval_nvmix_integral_c(double *lower, double *upper, double *U, int n, int d, int r,
                             int *kfactor, double *factor, double ZERO, double ONE, double *res)
{
    double yorg[r-1], sqrtmixorg, dorg, difforg, scprodorg;
    double yant[r-1], sqrtmixant, dant, diffant, scprodant;
    double ldorg, ldant, ldifforg, ldiffant, lforg, lfant;
    /* Note: <name>org stands for "original", <name>ant for antithetic */
    /* l<name>org/ant stands for log of that variable */
    /* y:       vector to save phi^{-1}(dj+uj(ej-dj)) */
    /* sqrtmix: used to store sqrt(F_w^{-1}(u_0)) */
    /* d:       current values of di from the paper */
    /* diff:    current values of (ei-di) from the paper */
    /* f:       current value of (e1-d1) * (e2-d2) * ... * (ei-di) */
    /* scprod:  scalar product sum factor_{ij} y_j */
    
    double tmp; /* to store temporary values */
    double lowermaxorg, upperminorg, lowermaxant, upperminant; /* needed in singular case */
    double scprodorgnew, scprodantnew;
    double sum = 0; /* to store sum_{i=1}^n (y_i + yant_i)/2 */
    double sumsq = 0; /* to store sum_{i=1}^n (y_i + yant_i)^2/4 (for variance calculation) */
    int current_limit; /* index of current element in 'lower'/'upper' */
    int i, j, l, m; /* counters for loops */
    
    /* Avoid recalculation */
    double qnormzero = qnorm(ZERO, 0, 1, 1, 0);
    double qnormone  = -qnormzero;

    /* For each row of U (for access, use U[s,k] = U[k * numrows + s]) */
    for(j = 0; j < n; j++){
        /* initialize current_limit */
        current_limit = 0;
        /* Grab the realizations of sqrt(mix) = sqrt(W) */
        sqrtmixorg = U[j];
        sqrtmixant = U[r*n+j];
        
        /* Grab 'lower' and 'upper' */
        lowermaxorg = lower[current_limit];
        upperminorg = upper[current_limit];
    
        /* Non-singular case: */
        if(r == d){
            lowermaxorg = lowermaxorg / factor[0];
            upperminorg = upperminorg / factor[0];
        } else if(kfactor[0] > 1){
            /* Singular case: Find active limit */
            for(i = 1; i <= kfactor[0]; i++){
                if(lower[i] > lowermaxorg){
                    lowermaxorg = lower[i];
                }
                if(upper[i] < upperminorg){
                    upperminorg = upper[i];
                }
            }
        }
        current_limit += kfactor[0];
        
        
        /* Need essentially log( pnorm(upper*) - pnorm(lower*)) */
        /* For higher accuracy, make use of 'lower.tail' argument in pnorm() when helpful */
        /* Note the arguments: .Call(C_pnorm, q, mean, sd, lower.tail, log.p) */
        
        /* lower_ = -Inf */
        if(lowermaxorg == R_NegInf){
            dorg = 0;
            dant = 0;
            if(upperminorg == R_PosInf){
                /* Case Phi(Inf) - Phi(-Inf) */
                difforg = 1;
                diffant = 1;
                ldifforg = 0;
                ldiffant = 0;
            } else {
                /* Case Phi(upper) - Phi(-Inf) = Phi(upper) */
                ldifforg = pnorm(upperminorg / sqrtmixorg, 0, 1, 1, 1);
                ldiffant = pnorm(upperminorg / sqrtmixant, 0, 1, 1, 1);
                difforg = exp(ldifforg);
                diffant = exp(ldiffant);
            }
        } else {
            /* lower_ != - Inf */
            if(upperminorg == R_PosInf){
                /* Case Phi(Inf) - Phi(lower_) = Phi(lower_, lower.tail = FALSE) */
                ldifforg = pnorm(lowermaxorg / sqrtmixorg, 0, 1, 0, 1);
                ldiffant = pnorm(lowermaxorg / sqrtmixant, 0, 1, 0, 1);
                difforg  = exp(ldifforg);
                diffant  = exp(ldiffant);
                dorg = 1 - difforg;
                dant = 1 - diffant;
            } else {
                /* Case Phi(upper_) - Phi(lower_)  */
                ldorg = pnorm(lowermaxorg / sqrtmixorg, 0, 1, 1, 1);
                ldant = pnorm(lowermaxorg / sqrtmixant, 0, 1, 1, 1);
                dorg  = exp(ldorg);
                dant  = exp(ldant);
                /* logsumexp trick for log(Phi(upper_)-Phi(lower_))*/
                tmp = pnorm(upperminorg / sqrtmixorg, 0, 1, 1, 1);
                ldifforg = tmp + log1p(-exp(ldorg - tmp));
                tmp = pnorm(upperminorg / sqrtmixant, 0, 1, 1, 1);
                ldiffant = tmp + log1p(-exp(ldant - tmp));
                difforg = exp(ldifforg);
                diffant = exp(ldiffant);
            }
        }

        /* Go through all r-1 columns (without first and last) */
        /* For better readability, we start at i = 0 */
        lforg = ldifforg;
        lfant = ldiffant;
        for(i = 0; i < r-1; i++){
            /* U[i * n + j] corresponds to U[j,i] in the orginal matrix */
            tmp = dorg + U[(i+1) * n + j] * difforg;
            /* Check if too close to 0 or 1 */
            if(tmp < ZERO){
                yorg[i] = qnormzero;
            } else if(tmp > ONE){
                yorg[i] = qnormone;
            } else {
                yorg[i] = qnorm(tmp, 0, 1, 1, 0);
            }

            /* The same for the antithetic value */
            tmp = dant + (1-U[(i+1) * n + j]) * diffant;
            if(tmp < ZERO){
                yant[i] = qnormzero;
            } else if(tmp > ONE){
                yant[i] = qnormone;
            } else {
                yant[i] = qnorm(tmp, 0, 1, 1, 0);
            }

            
            /* Calculate the scalar product sum factor[i,j] y[j] for j = 1:(i-1) */
            scprodorg = 0;
            scprodant = 0;
            for(l = 0; l < (i+1); l++){
                /* factor[l * d + current_limit] corresponds to factor[current_limit,l] in the original matrix */
                scprodorg += yorg[l] * factor[l * d + current_limit];
                scprodant += yant[l] * factor[l * d + current_limit];
            }
            
            lowermaxorg = (lower[current_limit] / sqrtmixorg - scprodorg);
            upperminorg = (upper[current_limit] / sqrtmixorg - scprodorg);
            lowermaxant = (lower[current_limit] / sqrtmixant - scprodant);
            upperminant = (upper[current_limit] / sqrtmixant - scprodant);
            
            /* Non-singular case: */
            if(r == d){
                /* Divide by C[i,i] != 0 */
                lowermaxorg = lowermaxorg / factor[current_limit * (d+1)];
                upperminorg = upperminorg / factor[current_limit * (d+1)];
                lowermaxant = lowermaxant / factor[current_limit * (d+1)];
                upperminant = upperminant / factor[current_limit * (d+1)];
            } else if(kfactor[i+1] > 1){
                /* Singular case: Go through the next rows of factor, adjust limit. Note: In the
                 singular case, 'factor's right-most element is always 1. */
                for(m = 1; m < (kfactor[i+1]+1); m++){
                    /* Calculate scalarproduct factor[i+l,j] y[j] for j = 1:(i-1) (next row of 'factor') */
                    scprodorgnew = 0;
                    scprodantnew = 0;
                    for(l = 0; l < (i+1); l++){
                        /* factor[l * d + i+1] corresponds to factor[i+1,l] in the original matrix */
                        scprodorgnew += yorg[l] * factor[l * d + current_limit + m];
                        scprodantnew += yant[l] * factor[l * d + current_limit + m];
                    }
                    /* Update 'lowermaxorg' if necessary */
                    tmp = (lower[current_limit + m] / sqrtmixorg - scprodorgnew);
                    if(tmp > lowermaxorg){
                        lowermaxorg = tmp;
                    }
                    /* Update 'upperminorg' if necessary */
                    tmp = (upper[current_limit + m] / sqrtmixorg - scprodorgnew);
                    if(tmp < upperminorg){
                        upperminorg = tmp;
                    }
                    /* Update 'lowermaxant' if necessary */
                    tmp = (lower[current_limit + m] / sqrtmixant - scprodantnew);
                    if(tmp > lowermaxant){
                        lowermaxant = tmp;
                    }
                    /* Update 'upperminant' if necessary */
                    tmp = (upper[current_limit + m] / sqrtmixant - scprodantnew);
                    if(tmp < upperminant){
                        upperminant = tmp;
                    }
                }
            }
            /* Calculate new d */
            /* Note: lower/upper<...>org Inf/-Inf <=> lower/upper<...>ant Inf/-Inf*/
            /* lower_ = -Inf */
            if(lowermaxorg == R_NegInf){
                dorg = 0;
                dant = 0;
                if(upperminorg == R_PosInf){
                    /* Case Phi(Inf) - Phi(-Inf) */
                    difforg = 1;
                    diffant = 1;
                    ldifforg = 0;
                    ldiffant = 0;
                } else {
                    /* Case Phi(upper) - Phi(-Inf) = Phi(upper) */
                    ldifforg = pnorm(upperminorg, 0, 1, 1, 1);
                    ldiffant = pnorm(upperminant, 0, 1, 1, 1);
                    difforg = exp(ldifforg);
                    diffant = exp(ldiffant);
                }
            } else {
                /* lower_ != -Inf */
                if(upperminorg == R_PosInf){
                    /* Case Phi(Inf) - Phi(lower_) = Phi(lower_, lower.tail = FALSE) */
                    ldifforg = pnorm(lowermaxorg, 0, 1, 0, 1);
                    ldiffant = pnorm(lowermaxant, 0, 1, 0, 1);
                    difforg  = exp(ldifforg);
                    diffant  = exp(ldiffant);
                    dorg = 1 - difforg;
                    dant = 1 - diffant;
                } else {
                    /* Case Phi(upper_) - Phi(lower_)  */
                    ldorg = pnorm(lowermaxorg, 0, 1, 1, 1);
                    ldant = pnorm(lowermaxant, 0, 1, 1, 1);
                    dorg  = exp(ldorg);
                    dant  = exp(ldant);
                    /* logsumexp trick for log(Phi(upper_)-Phi(lower_))*/
                    tmp = pnorm(upperminorg, 0, 1, 1, 1);
                    ldifforg = tmp + log1p(-exp(ldorg - tmp));
                    tmp = pnorm(upperminant, 0, 1, 1, 1);
                    ldiffant = tmp + log1p(-exp(ldant - tmp));
                    difforg = exp(ldifforg);
                    diffant = exp(ldiffant);
                }
                
            }
            /* Update products 'forg' and 'fant' */
            lforg += ldifforg;
            lfant += ldiffant;
            /* Update i in the singular case (as some rows may need to be skipped) */
            current_limit += kfactor[i+1];
        }
        tmp = (exp(lforg) + exp(lfant))/2;
        sum += tmp;
        sumsq += tmp*tmp;
    }
    res[0] = sum/n;
    res[1] = (sumsq - n * res[0] * res[0])/(n-1);
}

/**
 * @title R Interface for eval_nvmix_integral_c()
 * @param see eval_nvmix_integral_c()
 * @return mean(f(U)) where f is the integrand and U specifies the point-set
 * @author Erik Hintz, Marius Hofert (polishing)
 */
SEXP eval_nvmix_integral(SEXP lower, SEXP upper, SEXP U, SEXP n, SEXP d,
                         SEXP r, SEXP kfactor, SEXP factor, SEXP ZERO, SEXP ONE)
{
    SEXP res = PROTECT(allocVector(REALSXP, 2)); /* allocate memory */
    double *res_ = REAL(res); /* pointer to values of res */
    
    /* Main */
    eval_nvmix_integral_c(REAL(lower), REAL(upper), REAL(U), INTEGER(n)[0],
                          INTEGER(d)[0], INTEGER(r)[0], INTEGER(kfactor),
                          REAL(factor), REAL(ZERO)[0], REAL(ONE)[0], res_);
    UNPROTECT(1);
    /* Return */
    return res;
}
