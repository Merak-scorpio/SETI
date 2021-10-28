#pragma once

#include "constance.hpp"
#include "galaxy.hpp"


double IMF_rev(double phi);

void planet_loop(location *location);

void planet_nolife_process(location *location);

void planet_nointelligence_process(location *location);

void planet_intelligence_process(location *location);