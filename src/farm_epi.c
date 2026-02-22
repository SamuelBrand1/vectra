/**
 * @file farm_epi.c
 * @brief Per-farm epidemic steps: weather, deaths/recoveries,
 *        and bidirectional midge–host transmission.
 */

#include "simulation.h"
#include "simulation_internal.h"
#include "random.h"
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Forward declaration of control helper needed by detection logic */
static void implement_local_movement_ban(SimulationState *state, int centre_id,
                                         const ControlParams *ctrl);

/* ================================================================== */
/*  Get weather                                                        */
/* ================================================================== */

void farm_get_weather(Farm *farm, const SimulationState *state, gsl_rng *rng) {
    int day = state->simulation_day;
    farm->temp_today = state->temp_grid[farm->temp_grid_y][farm->temp_grid_x][day];
    farm->mean_rain_last_week = state->rain_grid[farm->rain_grid_y][farm->rain_grid_x][day];
    farm->wind_today = 0.0;
    farm->autocorr = 0.0;
    farm->overdispersion = (1.08 + 0.3763) * gsl_ran_gaussian(rng, 1.0);
}

/* ================================================================== */
/*  Deaths, recoveries, and detection                                  */
/* ================================================================== */

void farm_deaths_and_recoveries(Farm *farm, SimulationState *state, const EpiParams *epi,
                                const ControlParams *ctrl, gsl_rng *rng) {
    double dt_farm = 0.1; /* sub-day timestep */
    double sheep_mort_rate = 0.0055;

    /* Sheep recovery and mortality */
    if (num_inf_sheep(farm, epi->num_inf_stages_sheep) > 0.0) {
        double stages = (double)epi->num_inf_stages_sheep;
        for (double t = 0.0; t < 1.0; t += dt_farm) {
            /* Last stage: recovery */
            int last = epi->num_inf_stages_sheep - 1;
            int X = int_min(
                rand_poisson(rng, dt_farm * stages * epi->rec_rate_sheep * farm->i_sheep[last]),
                (int)farm->i_sheep[last]);
            farm->i_sheep[last] -= (double)X;
            farm->r_sheep += (double)X;

            /* Last stage: mortality */
            X = int_min(rand_poisson(rng, dt_farm * sheep_mort_rate * farm->i_sheep[last]),
                        (int)farm->i_sheep[last]);
            if (X > 0 && !farm->detected) {
                farm->detected = true;
                state->num_farms_detected_today++;
                if (!ctrl->no_control) {
                    if (!ctrl->no_farm_ban)
                        farm->movement_banned = true;
                    implement_local_movement_ban(state, farm->id, ctrl);
                    if (!state->btv_observed) {
                        state->btv_observed = true;
                        state->first_detected_farm_id = farm->id;
                    }
                }
            }
            farm->i_sheep[last] -= (double)X;
            state->num_sheep_deaths += X;

            /* Earlier stages: progression and mortality */
            for (int n = last - 1; n >= 0; n--) {
                X = int_min(
                    rand_poisson(rng, dt_farm * stages * epi->rec_rate_sheep * farm->i_sheep[n]),
                    (int)farm->i_sheep[n]);
                farm->i_sheep[n] -= (double)X;
                farm->i_sheep[n + 1] += (double)X;

                X = int_min(rand_poisson(rng, dt_farm * sheep_mort_rate * farm->i_sheep[n]),
                            (int)farm->i_sheep[n]);
                if (X > 0 && !farm->detected) {
                    farm->detected = true;
                    state->num_farms_detected_today++;
                    if (!ctrl->no_control) {
                        if (!ctrl->no_farm_ban)
                            farm->movement_banned = true;
                        implement_local_movement_ban(state, farm->id, ctrl);
                        if (!state->btv_observed) {
                            state->btv_observed = true;
                            state->first_detected_farm_id = farm->id;
                        }
                    }
                }
                farm->i_sheep[n] -= (double)X;
                state->num_sheep_deaths += X;
            }
        }
    }

    /* Cattle recovery (no mortality for cattle in this model) */
    if (num_inf_cattle(farm, epi->num_inf_stages_cattle) > 0.0) {
        double stages = (double)epi->num_inf_stages_cattle;
        for (double t = 0.0; t < 1.0; t += dt_farm) {
            int last = epi->num_inf_stages_cattle - 1;
            int X = int_min(
                rand_poisson(rng, dt_farm * stages * epi->rec_rate_cattle * farm->i_cattle[last]),
                (int)farm->i_cattle[last]);
            farm->i_cattle[last] -= (double)X;
            farm->r_cattle += (double)X;

            for (int n = last - 1; n >= 0; n--) {
                X = int_min(
                    rand_poisson(rng, dt_farm * stages * epi->rec_rate_cattle * farm->i_cattle[n]),
                    (int)farm->i_cattle[n]);
                farm->i_cattle[n] -= (double)X;
                farm->i_cattle[n + 1] += (double)X;
            }
        }
    }

    /* Passive detection */
    if (!farm->detected) {
        double inf_c = num_inf_cattle(farm, epi->num_inf_stages_cattle);
        double inf_s = num_inf_sheep(farm, epi->num_inf_stages_sheep);
        if (inf_c + inf_s > 0.0) {
            double not_detect_c = exp(inf_c * log(1.0 - epi->detection_prob_cattle));
            double not_detect_s = exp(inf_s * log(1.0 - epi->detection_prob_sheep));
            if (gsl_rng_uniform(rng) <= 1.0 - not_detect_c * not_detect_s) {
                farm->detected = true;
                state->num_farms_detected_today++;
                if (!ctrl->no_control) {
                    if (!ctrl->no_farm_ban)
                        farm->movement_banned = true;
                    implement_local_movement_ban(state, farm->id, ctrl);
                    if (!state->btv_observed) {
                        state->btv_observed = true;
                        state->first_detected_farm_id = farm->id;
                    }
                }
            }
        }
    }
}

/* ================================================================== */
/*  Midge-to-host transmission                                         */
/* ================================================================== */

void farm_transmission_midges_to_hosts(Farm *farm, SimulationState *state, const EpiParams *epi,
                                       const VectorSpecies *species, gsl_rng *rng) {
    double biting_rate = species->biting_rate(farm->temp_today);
    double biting_prob = 1.0 - exp(-biting_rate);
    double inf_density = state->inf_midge_density[farm->midge_grid_y][farm->midge_grid_x];
    double force = farm->rel_local_weight * inf_density * biting_prob;
    farm->force = force;

    double eff_animals = eff_num_animals(farm, epi->preference_for_sheep,
                                         epi->num_inf_stages_cattle, epi->num_inf_stages_sheep);
    if (eff_animals < 1.0)
        return;

    double prob_bite_sheep = epi->preference_for_sheep / eff_animals;
    double prob_bite_cattle = 1.0 / eff_animals;
    double prob_inf_sheep = 1.0 - exp(-force * prob_bite_sheep * epi->p_h);
    double prob_inf_cattle = 1.0 - exp(-force * prob_bite_cattle * epi->p_h);

    /* New sheep infections */
    int A;
    if (farm->s_sheep > 100.0 && prob_inf_sheep < 0.01 && farm->s_sheep * prob_inf_sheep < 20.0) {
        A = int_min(rand_poisson(rng, farm->s_sheep * prob_inf_sheep), (int)farm->s_sheep);
    } else {
        A = rand_binomial(rng, (int)farm->s_sheep, prob_inf_sheep);
    }

    /* New cattle infections */
    int B;
    if (farm->s_cattle > 100.0 && prob_inf_cattle < 0.01 &&
        farm->s_cattle * prob_inf_cattle < 20.0) {
        B = int_min(rand_poisson(rng, farm->s_cattle * prob_inf_cattle), (int)farm->s_cattle);
    } else {
        B = rand_binomial(rng, (int)farm->s_cattle, prob_inf_cattle);
    }

    farm->s_sheep -= (double)A;
    farm->i_sheep[0] += (double)A;
    state->num_sheep_infected_today += A;

    farm->s_cattle -= (double)B;
    farm->i_cattle[0] += (double)B;
    state->num_cattle_infected_today += B;
}

/* ================================================================== */
/*  Host-to-midge transmission                                         */
/* ================================================================== */

void farm_transmission_hosts_to_midges(const Farm *farm, SimulationState *state,
                                       const EpiParams *epi) {
    int day = state->simulation_day;

    /* Midge activity only during active season (day 60–330) */
    if (day % 365 <= 60 || day % 365 >= 330)
        return;

    double doy = (double)day;
    double climate = farm->v_intercept;
    climate += farm->sin_yearly * sin(2.0 * M_PI * doy / 365.25);
    climate += farm->cos_yearly * cos(2.0 * M_PI * doy / 365.25);
    climate += farm->sin_6_month * sin(4.0 * M_PI * doy / 365.25);
    climate += farm->cos_6_month * cos(4.0 * M_PI * doy / 365.25);
    climate += farm->cos_4_month * cos(6.0 * M_PI * doy / 365.25);
    climate +=
        farm->temp_eff * farm->temp_today + farm->temp_eff_sq * farm->temp_today * farm->temp_today;
    climate += farm->overdispersion + farm->autocorr;

    double bites_per_animal = epi->transmission_scalar * exp(climate);
    if (bites_per_animal > 5000.0)
        bites_per_animal = 5000.0;

    double eff_inf = eff_num_inf_animals(farm, epi->preference_for_sheep,
                                         epi->num_inf_stages_cattle, epi->num_inf_stages_sheep);

    double new_latent = epi->p_v * eff_inf * bites_per_animal;
    state->latent_midge_density[farm->midge_grid_y][farm->midge_grid_x][0] += new_latent;
}

/* ================================================================== */
/*  Local movement ban (needed by detection logic above)               */
/* ================================================================== */

static void implement_local_movement_ban(SimulationState *state, int centre_id,
                                         const ControlParams *ctrl) {
    Farm *centre = &state->farms[centre_id];

    /* Build local farm list on first detection */
    if (!centre->ever_been_detected) {
        centre->num_local_farms = 0;
        for (int k = 0; k < state->num_farms; k++) {
            if (k != centre_id &&
                dist_sq(&state->farms[k], centre) < ctrl->ban_radius * ctrl->ban_radius) {
                centre->local_farm_ids[centre->num_local_farms++] = k;
            }
        }
        centre->ever_been_detected = true;
    }

    if (!ctrl->no_farm_ban) {
        for (int n = 0; n < centre->num_local_farms; n++) {
            int bid = centre->local_farm_ids[n];
            state->farms[bid].movement_banned = true;
            state->farms[bid].free_area = false;
        }
    }

    /* County ban */
    if (ctrl->county_ban) {
        double cn = centre->county_number;
        for (int k = 0; k < state->num_farms; k++) {
            if (state->farms[k].county_number == cn) {
                state->farms[k].movement_banned = true;
                state->farms[k].free_area = false;
            }
        }
    }

    /* Total ban */
    if (ctrl->total_ban) {
        for (int k = 0; k < state->num_farms; k++) {
            state->farms[k].movement_banned = true;
            state->farms[k].free_area = false;
        }
    }
}
