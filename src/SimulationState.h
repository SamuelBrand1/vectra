/**
 * @file SimulationState.h
 * @brief Mutable simulation state that evolves during a model run.
 *
 * Contains the Farm struct and the top-level SimulationState that holds
 * all time-varying data: farm states, midge grids, daily counters, and
 * outbreak tracking flags.
 */

#ifndef SIMULATION_STATE_H
#define SIMULATION_STATE_H

#include <stdbool.h>

/** Maximum number of farms supported */
#define MAX_FARMS 200000

/** Maximum Erlang stages for infectious period */
#define MAX_INF_STAGES 20

/** Maximum stages for extrinsic incubation period */
#define MAX_EIP_STAGES 20

/** Maximum number of directed movement links across all farms */
#define MAX_MOVEMENT_LINKS 1000000

/** Grid dimensions (based on 5km cells over GB extent) */
#define MAX_GRID_S 244
#define MAX_GRID_E 131

/* ------------------------------------------------------------------ */
/*  Farm                                                               */
/* ------------------------------------------------------------------ */

/**
 * @brief State of a single farm during the simulation.
 */
typedef struct {
    int id;

    /* Location */
    double x[2];          /**< Coordinates in BNG (easting, northing) in metres */
    double county_number; /**< CPH county number */

    /* Grid cell indices for this farm */
    int temp_grid_x;
    int temp_grid_y;
    int rain_grid_x;
    int rain_grid_y;
    int midge_grid_x;
    int midge_grid_y;
    int ac_grid_x;
    int ac_grid_y;

    /* Midge abundance random effects (farm-level coefficients) */
    double v_intercept;
    double sin_yearly;
    double cos_yearly;
    double sin_6_month;
    double cos_6_month;
    double cos_4_month;
    double temp_eff;
    double temp_eff_sq;
    double rain_eff;
    double wind_eff;
    double autocorr;
    double overdispersion;

    /* Host demography and SIR state */
    double number_of_sheep;
    double number_of_cattle;
    double s_sheep;
    double i_sheep[MAX_INF_STAGES];
    double r_sheep;
    double s_cattle;
    double i_cattle[MAX_INF_STAGES];
    double r_cattle;

    /* Transmission */
    double rel_local_weight; /**< Relative attractiveness weight among local farms */
    double force;            /**< Force of infection on this farm */

    /* Control status */
    bool detected;
    bool movement_banned;
    bool protection_zone;
    bool surveillance_zone;
    bool free_area;
    bool ever_been_detected;
    bool ever_been_infected;
    bool first_infected_due_to_movement;

    /* Cached local farm list for movement bans */
    int num_local_farms;
    int local_farm_ids[MAX_FARMS];

    /* Today's weather at this farm */
    double temp_today;
    double mean_rain_last_week;
    double wind_today;
} Farm;

/* ------------------------------------------------------------------ */
/*  SimulationState                                                    */
/* ------------------------------------------------------------------ */

/**
 * @brief Top-level mutable state for a single simulation run.
 */
typedef struct {
    /* Time */
    int simulation_day; /**< Days elapsed since simulation start */
    int day_of_year;    /**< Current Julian day (1-365) */

    /* Farms */
    Farm farms[MAX_FARMS];
    int num_farms;

    /* Midge density grids */
    double latent_midge_density[MAX_GRID_S][MAX_GRID_E][MAX_EIP_STAGES];
    double inf_midge_density[MAX_GRID_S][MAX_GRID_E];
    double farm_biting_pref_grid[MAX_GRID_S][MAX_GRID_E];
    double diffusion_soln_grid[MAX_GRID_S][MAX_GRID_E];
    double diffusion_grid[MAX_GRID_S][MAX_GRID_E];

    /* Weather grids (read-only during simulation, loaded at setup) */
    double temp_grid[MAX_GRID_S][MAX_GRID_E][365];
    double rain_grid[MAX_GRID_S][MAX_GRID_E][365];
    double ac_grid[MAX_GRID_S][MAX_GRID_E];

    /* Movement network (sparse edge list) */
    int num_movement_links;                   /**< Total number of directed links */
    int movement_from[MAX_MOVEMENT_LINKS];    /**< Source farm ID for each link */
    int movement_to[MAX_MOVEMENT_LINKS];      /**< Destination farm ID for each link */
    double movement_risk[MAX_MOVEMENT_LINKS]; /**< Daily probability of movement for each link */

    /* Daily counters (reset each day) */
    int num_farms_detected_today;
    int num_sheep_infected_today;
    int num_cattle_infected_today;
    int num_sheep_deaths;

    /* Cumulative counters */
    int interrupted_movements;
    int days_of_movement_ban;
    int num_farms_checked;
    int num_tests;
    int num_pos_tests;
    int total_farm_days_movement_banned;
    int total_farm_days_affected_by_control;
    int num_movement_transmissions;
    int num_risky_moves_blocked;

    /* Outbreak tracking */
    bool btv_observed;
    int first_detected_farm_id;
    bool restriction_zones_implemented;
    bool active_surveillance_performed;
    int days_since_last_detection;
} SimulationState;

#endif /* SIMULATION_STATE_H */
