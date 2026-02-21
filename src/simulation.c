/**
 * @file simulation.c
 * @brief Top-level daily simulation orchestrator.
 *
 * Calls each simulation step in order. The step implementations
 * live in separate files:
 * - midge_dynamics.c — mortality, EIP incubation, diffusion
 * - movement.c       — livestock movement transmission
 * - control.c        — restriction zones, movement bans, surveillance
 * - farm_epi.c       — per-farm weather, SIR dynamics, transmission
 */

#include "simulation.h"

void simulate_day(SimulationState *state, const SimulationParams *sim, const EpiParams *epi,
                  const ControlParams *ctrl, const GridParams *grids, const VectorSpecies *species,
                  gsl_rng *rng) {
    /* Reset daily counters */
    state->num_farms_detected_today = 0;
    state->num_sheep_infected_today = 0;
    state->num_cattle_infected_today = 0;
    state->num_sheep_deaths = 0;

    /* 1. Control measures */
    apply_control_measures(state, epi, ctrl);

    /* 2. Midge mortality and EIP progression */
    midge_mortality_and_incubation(state, epi, grids, species);

    /* 3. Midge diffusion */
    midge_diffusion_for_day(state, sim, epi, grids);

    /* 4. Livestock movement transmission */
    movement_transmission(state, epi, ctrl, rng);

    /* 5. Per-farm epidemic updates */
    for (int k = 0; k < state->num_farms; k++) {
        farm_get_weather(&state->farms[k], state, rng);
        farm_deaths_and_recoveries(&state->farms[k], state, epi, ctrl, rng);
        farm_transmission_midges_to_hosts(&state->farms[k], state, epi, species, rng);
        farm_transmission_hosts_to_midges(&state->farms[k], state, epi);
    }

    /* 6. Advance time */
    state->simulation_day++;
    state->day_of_year = state->simulation_day % 365;
}
