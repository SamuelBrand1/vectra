/**
 * @file random.h
 * @brief GSL wrapper functions for random number generation.
 *
 * Thin wrappers around GSL random distribution functions used by the
 * simulation. All functions take an explicit gsl_rng pointer — no
 * global generator state.
 *
 * Why wrap GSL rather than calling it directly?
 * - Consistent snake_case naming across the codebase.
 * - Argument-order safety: GSL's binomial takes (rng, p, n) which is
 *   easy to confuse; our wrapper takes the more natural (rng, n, p).
 * - Single place to change if GSL is ever swapped for another library.
 * - Houses rand_neg_binomial(), which implements the Poisson–Gamma
 *   mixture rather than delegating to a single GSL call.
 */

#ifndef RANDOM_H
#define RANDOM_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/** @brief Draw from Binomial(n, p). */
int rand_binomial(const gsl_rng *rng, int n, double p);

/** @brief Draw from Poisson(lambda). */
int rand_poisson(const gsl_rng *rng, double lambda);

/** @brief Draw from Gamma(shape, scale). */
double rand_gamma(const gsl_rng *rng, double shape, double scale);

/** @brief Draw from NegativeBinomial(k, p) via Poisson–Gamma mixture. */
int rand_neg_binomial(const gsl_rng *rng, double k, double p);

/** @brief Poisson PMF: P(X = x) for X ~ Poisson(lambda). */
double poisson_pmf(int x, double lambda);

/** @brief Poisson CDF: P(X <= x) for X ~ Poisson(lambda). */
double poisson_cdf(int x, double lambda);

/** @brief Poisson survival: P(X > x) for X ~ Poisson(lambda). */
double poisson_sf(int x, double lambda);

/** @brief Return the smaller of a and b. */
int int_min(int a, int b);

#endif /* RANDOM_H */
