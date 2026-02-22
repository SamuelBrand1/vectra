/**
 * @file control.c
 * @brief Disease control measures: movement bans, restriction zones,
 *        and active surveillance.
 */

#include "simulation.h"
#include "simulation_internal.h"

/**
 * @brief Find farms within ban_radius of centre and ban their movements.
 */
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

/**
 * @brief Set up protection and surveillance zones around a farm.
 */
static void setup_restriction_zone(SimulationState *state, int centre_id,
                                   const ControlParams *ctrl) {
    Farm *centre = &state->farms[centre_id];
    double pz2 = ctrl->pz_radius * ctrl->pz_radius;
    double sz2 = ctrl->sz_radius * ctrl->sz_radius;

    for (int k = 0; k < state->num_farms; k++) {
        double d2 = dist_sq(&state->farms[k], centre);
        if (d2 <= pz2) {
            state->farms[k].protection_zone = true;
            state->farms[k].free_area = false;
        } else if (d2 <= sz2) {
            state->farms[k].surveillance_zone = true;
            state->farms[k].free_area = false;
        }
    }
    state->restriction_zones_implemented = true;
}

/**
 * @brief Perform active surveillance around first detected farm.
 */
static void perform_active_surveillance(SimulationState *state, const EpiParams *epi) {
    double surv_radius = 15000.0;
    double sr2 = surv_radius * surv_radius;
    Farm *centre = &state->farms[state->first_detected_farm_id];

    for (int k = 0; k < state->num_farms; k++) {
        if (dist_sq(&state->farms[k], centre) <= sr2) {
            Farm *f = &state->farms[k];
            state->num_farms_checked++;
            state->num_tests += (int)(num_cattle(f, epi->num_inf_stages_cattle) +
                                      num_sheep(f, epi->num_inf_stages_sheep));
            double inf = num_inf_cattle(f, epi->num_inf_stages_cattle) +
                         num_inf_sheep(f, epi->num_inf_stages_sheep);
            if (inf > 0.0) {
                f->detected = true;
                state->num_pos_tests += (int)(inf + f->r_sheep + f->r_cattle);
            }
        }
    }
    state->active_surveillance_performed = true;
}

void apply_control_measures(SimulationState *state, const EpiParams *epi,
                            const ControlParams *ctrl) {
    if (ctrl->no_control)
        return;

    if (state->btv_observed && !state->restriction_zones_implemented) {
        if (ctrl->restriction_zones) {
            setup_restriction_zone(state, state->first_detected_farm_id, ctrl);
        }
        if (!state->active_surveillance_performed) {
            perform_active_surveillance(state, epi);
        }
    }
}
