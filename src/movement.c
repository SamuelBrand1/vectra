/**
 * @file movement.c
 * @brief Livestock movement transmission between farms.
 */

#include "simulation.h"
#include "simulation_internal.h"
#include "random.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

/**
 * @brief Process a single movement link between two farms.
 */
static void transmission_via_movement(SimulationState *state, int from_id, int to_id, double risk,
                                      const EpiParams *epi, const MovementParams *mov,
                                      const ControlParams *ctrl, gsl_rng *rng) {
    Farm *src = &state->farms[from_id];
    Farm *dst = &state->farms[to_id];

    /* Check if movement occurs today */
    double u = gsl_rng_uniform(rng);
    if (u > risk)
        return;

    /* Check for control-based interruption */
    bool interrupt = false;
    if (src->movement_banned || dst->movement_banned)
        interrupt = true;
    if (src->protection_zone && !dst->protection_zone)
        interrupt = true;
    if (src->surveillance_zone && dst->free_area)
        interrupt = true;

    if (interrupt) {
        state->interrupted_movements++;
        if (num_inf_cattle(src, epi->num_inf_stages_cattle) +
                num_inf_sheep(src, epi->num_inf_stages_sheep) >
            0) {
            state->num_risky_moves_blocked++;
        }
        return;
    }

    /* Determine cattle vs sheep move */
    int num_inf_moved = 0;
    double total_sheep = num_sheep(src, epi->num_inf_stages_sheep);
    double total_cattle = num_cattle(src, epi->num_inf_stages_cattle);

    if (total_sheep + total_cattle < 1.0)
        return;

    bool cattle_move = gsl_rng_uniform(rng) > total_sheep / (total_sheep + total_cattle);

    if (cattle_move) {
        int size =
            int_min(1 + (int)gsl_ran_negative_binomial(rng, mov->cattle_shipment_size_k,
                                                        mov->cattle_shipment_size_p),
                    (int)total_cattle);
        double inf_cattle = num_inf_cattle(src, epi->num_inf_stages_cattle);
        for (int n = 0; n < size; n++) {
            double density = inf_cattle / total_cattle;
            if (gsl_rng_uniform(rng) < density) {
                double sel = gsl_rng_uniform(rng) * inf_cattle;
                double cum = 0.0;
                for (int k = 0; k < epi->num_inf_stages_cattle; k++) {
                    cum += src->i_cattle[k];
                    if (cum >= sel) {
                        src->i_cattle[k] -= 1.0;
                        dst->i_cattle[k] += 1.0;
                        num_inf_moved++;
                        inf_cattle -= 1.0;
                        break;
                    }
                }
            }
        }
    } else {
        int size =
            int_min(1 + (int)gsl_ran_negative_binomial(rng, mov->sheep_shipment_size_k,
                                                        mov->sheep_shipment_size_p),
                    (int)total_sheep);
        double inf_sheep = num_inf_sheep(src, epi->num_inf_stages_sheep);
        for (int n = 0; n < size; n++) {
            double density = inf_sheep / total_sheep;
            if (gsl_rng_uniform(rng) < density) {
                double sel = gsl_rng_uniform(rng) * inf_sheep;
                double cum = 0.0;
                for (int k = 0; k < epi->num_inf_stages_sheep; k++) {
                    cum += src->i_sheep[k];
                    if (cum >= sel) {
                        src->i_sheep[k] -= 1.0;
                        dst->i_sheep[k] += 1.0;
                        num_inf_moved++;
                        inf_sheep -= 1.0;
                        break;
                    }
                }
            }
        }
    }

    if (num_inf_moved > 0) {
        state->num_movement_transmissions++;
    }
}

void movement_transmission(SimulationState *state, const EpiParams *epi, const MovementParams *mov,
                           const ControlParams *ctrl, gsl_rng *rng) {
    for (int k = 0; k < state->num_movement_links; k++) {
        transmission_via_movement(state, state->movement_from[k], state->movement_to[k],
                                  state->movement_risk[k], epi, mov, ctrl, rng);
    }
}
