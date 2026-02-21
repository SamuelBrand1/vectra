/**
 * @file entomology.c
 * @brief Implementations of vector species temperature-dependent rate functions.
 */

#include "entomology.h"
#include <math.h>

double culicoides_biting_rate(double temperature) {
    if (temperature > 3.7 && temperature < 41.9) {
        return 0.0002 * temperature * (temperature - 3.7) * exp(0.37 * log(41.9 - temperature));
    }
    return 0.0;
}

double culicoides_mortality_rate(double temperature) {
    if (temperature > -2.0) {
        return 0.009 * exp(0.16 * temperature);
    }
    return 100.0;
}

double culicoides_incubation_rate(double temperature) {
    double rate = 0.018 * (temperature - 13.4);
    if (rate > 0.0) {
        return rate;
    }
    return 0.0;
}

VectorSpecies culicoides_species(void) {
    VectorSpecies species = {
        .name = "Culicoides",
        .biting_rate = culicoides_biting_rate,
        .mortality_rate = culicoides_mortality_rate,
        .incubation_rate = culicoides_incubation_rate,
    };
    return species;
}
