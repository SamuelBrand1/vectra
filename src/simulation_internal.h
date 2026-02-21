/**
 * @file simulation_internal.h
 * @brief Shared helper functions used across simulation .c files.
 *
 * Not part of the public API — only included by the simulation
 * implementation files (control.c, farm_epi.c, movement.c, etc.).
 * All functions are static inline so each translation unit gets
 * its own copy with no linker conflicts.
 */

#ifndef SIMULATION_INTERNAL_H
#define SIMULATION_INTERNAL_H

#include "SimulationState.h"
#include "parameters.h"

static inline double dist_sq(const Farm *a, const Farm *b) {
    double dx = a->x[0] - b->x[0];
    double dy = a->x[1] - b->x[1];
    return dx * dx + dy * dy;
}

static inline double num_inf_cattle(const Farm *f, int num_stages) {
    double sum = 0.0;
    for (int i = 0; i < num_stages; i++) {
        sum += f->i_cattle[i];
    }
    return sum;
}

static inline double num_inf_sheep(const Farm *f, int num_stages) {
    double sum = 0.0;
    for (int i = 0; i < num_stages; i++) {
        sum += f->i_sheep[i];
    }
    return sum;
}

static inline double num_cattle(const Farm *f, int num_stages) {
    return f->s_cattle + num_inf_cattle(f, num_stages) + f->r_cattle;
}

static inline double num_sheep(const Farm *f, int num_stages) {
    return f->s_sheep + num_inf_sheep(f, num_stages) + f->r_sheep;
}

static inline double eff_num_animals(const Farm *f, double pref, int inf_stages_cattle,
                                     int inf_stages_sheep) {
    return num_cattle(f, inf_stages_cattle) + pref * num_sheep(f, inf_stages_sheep);
}

static inline double eff_num_inf_animals(const Farm *f, double pref, int inf_stages_cattle,
                                         int inf_stages_sheep) {
    double sum = 0.0;
    for (int i = 0; i < inf_stages_cattle; i++) {
        sum += f->i_cattle[i];
    }
    for (int i = 0; i < inf_stages_sheep; i++) {
        sum += pref * f->i_sheep[i];
    }
    return sum;
}

#endif /* SIMULATION_INTERNAL_H */
