/**
 * @file midge_dynamics.c
 * @brief Midge mortality, EIP incubation, and spatial diffusion.
 */

#include "simulation.h"
#include "random.h"
#include <math.h>

/* ================================================================== */
/*  Midge mortality and EIP incubation                                 */
/* ================================================================== */

void midge_mortality_and_incubation(SimulationState *state, const EpiParams *epi,
                                    const GridParams *grids, const VectorSpecies *species) {
    int R = (int)(grids->midge_grid_width / grids->temp_grid_width);
    int num_eip = epi->num_eip_stages;
    double soln[MAX_EIP_STAGES + 1];

    for (int i = 0; i < MAX_GRID_S; i += R) {
        for (int j = 0; j < MAX_GRID_E; j += R) {
            double temp = state->temp_grid[i][j][state->simulation_day];
            double mort = exp(-species->mortality_rate(temp));
            double incub = (double)num_eip * species->incubation_rate(temp);

            /* Apply mortality */
            state->inf_midge_density[i][j] *= mort;
            double sum = 0.0;
            for (int k = 0; k < num_eip; k++) {
                state->latent_midge_density[i][j][k] *= mort;
                sum += state->latent_midge_density[i][j][k];
            }

            /* EIP stage progression (Poisson-distributed transitions) */
            if (incub > 0.0 && sum > 0.0) {
                for (int n = 0; n < num_eip; n++) {
                    soln[n] = 0.0;
                    for (int k = 0; k <= n; k++) {
                        soln[n] += state->latent_midge_density[i][j][k] * poisson_pmf(n - k, incub);
                    }
                }
                /* Transitions to infectious state */
                soln[num_eip] = state->inf_midge_density[i][j];
                for (int k = 0; k < num_eip; k++) {
                    soln[num_eip] +=
                        state->latent_midge_density[i][j][k] * poisson_sf(num_eip - k - 1, incub);
                }
                /* Write back */
                for (int n = 0; n < num_eip; n++) {
                    state->latent_midge_density[i][j][n] = soln[n];
                }
                state->inf_midge_density[i][j] = soln[num_eip];
            }
        }
    }
}

/* ================================================================== */
/*  Midge diffusion                                                    */
/* ================================================================== */

/**
 * @brief Single diffusion sub-step for all midge grids.
 */
static void midge_diffusion_step(SimulationState *state, const SimulationParams *sim,
                                 const EpiParams *epi, const GridParams *grids) {
    double dt = sim->dt;
    double h2 = grids->midge_grid_width * grids->midge_grid_width;
    int num_eip = epi->num_eip_stages;

    /* Latent midge stages */
    for (int k = 0; k < num_eip; k++) {
        for (int i = 1; i < MAX_GRID_S - 1; i++) {
            for (int j = 1; j < MAX_GRID_E - 1; j++) {
                if (state->latent_midge_density[i][j][k] > 1e-5) {
                    double D = state->diffusion_grid[i][j];
                    double flux = D * dt * state->latent_midge_density[i][j][k] / h2;
                    state->diffusion_soln_grid[i][j] -= 2.0 * flux;
                    state->diffusion_soln_grid[i + 1][j] += 0.5 * flux;
                    state->diffusion_soln_grid[i - 1][j] += 0.5 * flux;
                    state->diffusion_soln_grid[i][j + 1] += 0.5 * flux;
                    state->diffusion_soln_grid[i][j - 1] += 0.5 * flux;
                }
            }
        }
        for (int i = 1; i < MAX_GRID_S - 1; i++) {
            for (int j = 1; j < MAX_GRID_E - 1; j++) {
                state->latent_midge_density[i][j][k] += state->diffusion_soln_grid[i][j];
                state->diffusion_soln_grid[i][j] = 0.0;
            }
        }
    }

    /* Infectious midges */
    for (int i = 1; i < MAX_GRID_S - 1; i++) {
        for (int j = 1; j < MAX_GRID_E - 1; j++) {
            if (state->inf_midge_density[i][j] > 1e-5) {
                double D = state->diffusion_grid[i][j];
                double flux = D * dt * state->inf_midge_density[i][j] / h2;
                state->diffusion_soln_grid[i][j] -= 2.0 * flux;
                state->diffusion_soln_grid[i + 1][j] += 0.5 * flux;
                state->diffusion_soln_grid[i - 1][j] += 0.5 * flux;
                state->diffusion_soln_grid[i][j + 1] += 0.5 * flux;
                state->diffusion_soln_grid[i][j - 1] += 0.5 * flux;
            }
        }
    }
    for (int i = 1; i < MAX_GRID_S - 1; i++) {
        for (int j = 1; j < MAX_GRID_E - 1; j++) {
            state->inf_midge_density[i][j] += state->diffusion_soln_grid[i][j];
            state->diffusion_soln_grid[i][j] = 0.0;
        }
    }
}

void midge_diffusion_for_day(SimulationState *state, const SimulationParams *sim,
                             const EpiParams *epi, const GridParams *grids) {
    double elapsed = 0.0;
    while (elapsed < 1.0) {
        midge_diffusion_step(state, sim, epi, grids);
        elapsed += sim->dt;
    }
}
