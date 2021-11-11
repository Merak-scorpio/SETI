#pragma once

#include "constance.hpp"
#include "galaxy.hpp"
#include "planet.hpp"

void SNe_loop(vector<SNe *> *SNe_list);

bool SNe_add_judge(location *location, SNe *SNe);

void SNe_add(location *location, vector<SNe *> *SNe_list);

bool excute_judge(double x1, double y1, double z1, double x2, double y2, double z2, double r);

// void SNe_excute(location *location, int SNe_num);
void SNe_excute(location *location);