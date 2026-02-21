/**
 * @file simulation.h
 * @brief Simulation step function declarations for the VECTRA model.
 *
 * Each function implements one step of the daily simulation loop. Functions
 * take a mutable SimulationState pointer and const parameter structs — no
 * global state is accessed. Stochastic functions take an explicit gsl_rng*.
 *
 * The daily simulation order is:
 * 1. apply_control_measures()
 * 2. midge_mortality_and_incubation()
 * 3. midge_diffusion_for_day()
 * 4. movement_transmission()
 * 5. Per farm: farm_get_weather(), farm_deaths_and_recoveries(),
 *    farm_transmission_midges_to_hosts(), farm_transmission_hosts_to_midges()
 * 6. Increment simulation day
 */

#ifndef SIMULATION_H
#define SIMULATION_H

#include "SimulationState.h"
#include "parameters.h"
#include "entomology.h"
#include <gsl/gsl_rng.h>

/* ------------------------------------------------------------------ */
/*  Top-level orchestration                                            */
/* ------------------------------------------------------------------ */

/**
 * @brief Run one complete day of the simulation.
 *
 * Calls all daily steps in order: control measures, midge dynamics,
 * movement transmission, and per-farm epidemic updates. Increments
 * the simulation day counter at the end.
 *
 * @param state   Mutable simulation state
 * @param sim     Simulation run parameters (timesteps, etc.)
 * @param epi     Epidemiological parameters
 * @param ctrl    Control measure parameters
 * @param grids   Spatial grid parameters
 * @param species Vector species (temperature-dependent rate functions)
 * @param rng     GSL random number generator
 */
void simulate_day(SimulationState *state, const SimulationParams *sim, const EpiParams *epi,
                  const ControlParams *ctrl, const GridParams *grids, const VectorSpecies *species,
                  gsl_rng *rng);

/* ------------------------------------------------------------------ */
/*  Grid-level steps                                                   */
/* ------------------------------------------------------------------ */

/**
 * @brief Apply temperature-dependent mortality and EIP progression to midges.
 *
 * For each grid cell, reduces latent and infectious midge densities by the
 * local mortality rate, then progresses latent midges through EIP stages
 * using Poisson-distributed stage transitions.
 *
 * @param state   Mutable simulation state (midge grids and temp grid)
 * @param epi     Epidemiological parameters (num_eip_stages)
 * @param grids   Grid parameters (grid widths for temp-to-midge mapping)
 * @param species Vector species (mortality and incubation rate functions)
 */
void midge_mortality_and_incubation(SimulationState *state, const EpiParams *epi,
                                    const GridParams *grids, const VectorSpecies *species);

/**
 * @brief Diffuse midge populations spatially for one day.
 *
 * Runs the 2D forward-in-time diffusion scheme in sub-steps of size dt
 * until one full day has elapsed. Applies to both latent (per EIP stage)
 * and infectious midge density grids.
 *
 * @param state Mutable simulation state (midge and diffusion grids)
 * @param sim   Simulation parameters (dt timestep)
 * @param epi   Epidemiological parameters (num_eip_stages)
 * @param grids Grid parameters (midge_grid_width)
 */
void midge_diffusion_for_day(SimulationState *state, const SimulationParams *sim,
                             const EpiParams *epi, const GridParams *grids);

/* ------------------------------------------------------------------ */
/*  System-level steps                                                 */
/* ------------------------------------------------------------------ */

/**
 * @brief Process all livestock movements and associated disease transmission.
 *
 * Iterates over the movement edge list. For each link, stochastically
 * determines whether a movement occurs, checks for control-based
 * interruption, and transfers infected animals if the movement proceeds.
 *
 * @param state Mutable simulation state (farms, movement network, counters)
 * @param epi   Epidemiological parameters (num_inf_stages)
 * @param ctrl  Control parameters (for movement ban checks)
 * @param rng   GSL random number generator
 */
void movement_transmission(SimulationState *state, const EpiParams *epi, const ControlParams *ctrl,
                           gsl_rng *rng);

/**
 * @brief Check outbreak detection flags and apply control measures.
 *
 * If BTV has been observed but restriction zones have not yet been set up,
 * establishes protection and surveillance zones around the first detected
 * farm. Performs active surveillance if applicable.
 *
 * @param state Mutable simulation state (farms, outbreak flags)
 * @param epi   Epidemiological parameters (for active surveillance counts)
 * @param ctrl  Control parameters (zone radii, ban policies)
 */
void apply_control_measures(SimulationState *state, const EpiParams *epi,
                            const ControlParams *ctrl);

/* ------------------------------------------------------------------ */
/*  Per-farm steps                                                     */
/* ------------------------------------------------------------------ */

/**
 * @brief Load today's weather data for a farm from the spatial grids.
 *
 * Looks up temperature and rainfall from the grid cell corresponding to
 * this farm's location and stores them on the farm struct.
 *
 * @param farm  The farm to update
 * @param state Simulation state (weather grids, simulation_day)
 * @param rng   GSL random number generator (for overdispersion sampling)
 */
void farm_get_weather(Farm *farm, const SimulationState *state, gsl_rng *rng);

/**
 * @brief Process animal deaths, recoveries, and disease detection on a farm.
 *
 * Progresses infected animals through Erlang stages (recovery), applies
 * sheep mortality, and checks for passive detection. If detected, triggers
 * movement bans and updates outbreak tracking flags on the state.
 *
 * @param farm  The farm to update
 * @param state Mutable simulation state (for counters and control triggers)
 * @param epi   Epidemiological parameters (recovery/mortality/detection rates)
 * @param ctrl  Control parameters (ban policies)
 * @param rng   GSL random number generator
 */
void farm_deaths_and_recoveries(Farm *farm, SimulationState *state, const EpiParams *epi,
                                const ControlParams *ctrl, gsl_rng *rng);

/**
 * @brief Transmit infection from infectious midges to susceptible livestock.
 *
 * Calculates the force of infection from the local infectious midge density
 * and the temperature-dependent biting rate, then stochastically infects
 * susceptible cattle and sheep.
 *
 * @param farm    The farm to update
 * @param state   Mutable simulation state (midge grids, daily counters)
 * @param epi     Epidemiological parameters (p_h, preference_for_sheep)
 * @param species Vector species (biting rate function)
 * @param rng     GSL random number generator
 */
void farm_transmission_midges_to_hosts(Farm *farm, SimulationState *state, const EpiParams *epi,
                                       const VectorSpecies *species, gsl_rng *rng);

/**
 * @brief Transmit infection from infected livestock to susceptible midges.
 *
 * Calculates the expected number of newly inoculated midges based on the
 * climate-driven midge abundance model and the number of infected animals,
 * then adds them to the first latent EIP stage on the local grid cell.
 *
 * @param farm  The farm to read from
 * @param state Mutable simulation state (latent midge grid)
 * @param epi   Epidemiological parameters (p_v, transmission_scalar, preference)
 */
void farm_transmission_hosts_to_midges(const Farm *farm, SimulationState *state,
                                       const EpiParams *epi);

#endif /* SIMULATION_H */
