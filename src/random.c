/**
 * @file random.c
 * @brief Implementations of GSL wrapper functions.
 */

#include "random.h"

int rand_binomial(const gsl_rng *rng, int n, double p) {
    return (int)gsl_ran_binomial(rng, p, (unsigned int)n);
}

int rand_poisson(const gsl_rng *rng, double lambda) {
    return (int)gsl_ran_poisson(rng, lambda);
}

double rand_gamma(const gsl_rng *rng, double shape, double scale) {
    return gsl_ran_gamma(rng, shape, scale);
}

int rand_neg_binomial(const gsl_rng *rng, double k, double p) {
    double g = rand_gamma(rng, k, p / (1.0 - p));
    return rand_poisson(rng, g);
}

double poisson_pmf(int x, double lambda) {
    return gsl_ran_poisson_pdf((unsigned int)x, lambda);
}

double poisson_cdf(int x, double lambda) {
    return gsl_cdf_poisson_P((unsigned int)x, lambda);
}

double poisson_sf(int x, double lambda) {
    return gsl_cdf_poisson_Q((unsigned int)x, lambda);
}

int int_min(int a, int b) {
    return (a < b) ? a : b;
}
