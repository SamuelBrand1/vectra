/**
 * @file parameters.h
 * @brief Parameter structs for the VECTRA BTV simulation model.
 *
 * Defines all configurable parameters grouped by function: simulation settings,
 * epidemiological rates, control measures, initial conditions, and spatial grids.
 */

#ifndef PARAMETERS_H
#define PARAMETERS_H

#include <stdbool.h>

/**
 * @brief Parameters controlling the simulation run.
 */
typedef struct {
    double dt;                 /**< Timestep for midge dynamics (days) */
    double dt_farm;            /**< Timestep for farm-level dynamics (days) */
    double initial_density_inf_midges;     /**< Initial density of infectious midges at seed location */
    double initial_width_from_central_site; /**< Spatial extent of initial seeding (metres) */
    double outbreak_county;                /**< County number where the outbreak is seeded */
    int num_days;              /**< Total number of days to simulate */
    int num_reps;              /**< Number of Monte Carlo repetitions */
    int start_day_of_year;     /**< Julian day to start (1-365) */
} SimulationParams;

/**
 * @brief Epidemiological parameters for BTV transmission and disease.
 */
typedef struct {
    double detection_prob_cattle;  /**< Daily probability of detecting an infected cow */
    double detection_prob_sheep;   /**< Daily probability of detecting an infected sheep */
    double diffusion_length_scale; /**< Length scale for midge diffusion (metres) */
    int num_inf_stages_sheep;      /**< Erlang stages for sheep infectious period */
    int num_inf_stages_cattle;     /**< Erlang stages for cattle infectious period */
    int num_eip_stages;            /**< Stages for extrinsic incubation period */
    double p_v;                    /**< Probability of vector infection per bite on infectious host */
    double p_h;                    /**< Probability of host infection per bite from infectious vector */
    double sheep_mort_rate;        /**< Daily mortality rate for infected sheep */
    double rec_rate_sheep;         /**< Recovery rate for sheep (per day) */
    double rec_rate_cattle;        /**< Recovery rate for cattle (per day) */
} EpiParams;

/**
 * @brief Parameters defining disease control measures.
 */
typedef struct {
    double ban_radius;       /**< Radius of local movement ban around detected farms (metres) */
    bool county_ban;         /**< Ban all movement within the county of a detected farm */
    bool no_control;         /**< Disable all control measures */
    bool no_farm_ban;        /**< Disable farm-level movement bans */
    bool pre_movement_tests; /**< Require pre-movement testing */
    double pz_radius;        /**< Radius of the protection zone (metres) */
    bool restriction_zones;  /**< Enable protection and surveillance zones */
    double sz_radius;        /**< Radius of the surveillance zone (metres) */
    bool total_ban;          /**< Ban all animal movement nationally */
} ControlParams;

/**
 * @brief Parameters defining the spatial grid resolutions.
 */
typedef struct {
    double autocorr_grid_width; /**< Grid cell width for autocorrelation field (metres) */
    double discretisation;      /**< Discretisation parameter for diffusion solver */
    double midge_grid_width;    /**< Grid cell width for midge density (metres) */
    double rain_grid_width;     /**< Grid cell width for rainfall data (metres) */
    double temp_grid_width;     /**< Grid cell width for temperature data (metres) */
} GridParams;

#endif /* PARAMETERS_H */