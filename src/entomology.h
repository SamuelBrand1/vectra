/**
 * @file entomology.h
 * @brief Vector species definitions and temperature-dependent rate functions.
 *
 * Uses function pointers to allow different vector species (e.g., Culicoides,
 * mosquitoes) to have different temperature-response curves. The simulation
 * code calls rates through the VectorSpecies struct, so adding a new species
 * only requires implementing three functions and wiring them up.
 */

#ifndef ENTOMOLOGY_H
#define ENTOMOLOGY_H

/** Function signature for temperature-dependent rate functions */
typedef double (*RateFunction)(double temperature);

/**
 * @brief Defines a vector species via its temperature-dependent rate functions.
 */
typedef struct {
    const char *name;              /**< Species name for logging */
    RateFunction biting_rate;      /**< Biting rate per day as f(T) */
    RateFunction mortality_rate;   /**< Daily mortality rate as f(T) */
    RateFunction incubation_rate;  /**< EIP progression rate as f(T) */
} VectorSpecies;

/* ------------------------------------------------------------------ */
/*  Culicoides (BTV midges) — default species                          */
/* ------------------------------------------------------------------ */

/** @brief Culicoides biting rate (per day) as a function of temperature. */
double culicoides_biting_rate(double temperature);

/** @brief Culicoides daily mortality rate as a function of temperature. */
double culicoides_mortality_rate(double temperature);

/** @brief Culicoides EIP progression rate as a function of temperature. */
double culicoides_incubation_rate(double temperature);

/**
 * @brief Returns a VectorSpecies configured for Culicoides.
 */
VectorSpecies culicoides_species(void);

#endif /* ENTOMOLOGY_H */
