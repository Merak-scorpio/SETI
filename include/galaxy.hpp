#pragma once

#include "constance.hpp"
#ifdef WIN32
#include <io.h>
#elif linux
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#endif

using namespace std;

struct nolife_planet {
  int    star_age      = 0;
  int    life_evo_time = 0;
  double x_coordinate;
  double y_coordinate;
  double z_coordinate;
  double star_T_L;
  double planet_T_min;
};

struct nointelligence_planet {
  int    star_age = 0;
  int    life_age = 0;
  double x_coordinate;
  double y_coordinate;
  double z_coordinate;
  double star_T_L;
};

struct intelligence_planet {
  int    star_age         = 0;
  int    intelligence_age = 0;
  double x_coordinate;
  double y_coordinate;
  double z_coordinate;
  double star_T_L;
};

struct SNe {
  double x_coordinate;
  double y_coordinate;
  double z_coordinate;
  double D_SNe;
};

struct location {
  int    x        = 500;
  int    y        = 500;
  double r        = 0;
  double Normal_A = 0;

  double fin_t         = 0;
  double sigma_SFR     = 0;
  double sigma_Gas     = 0;
  double Gas_mass_left = 0;
  double SFR_left      = 0;

  vector<nolife_planet *>         *Nolife_planets         = new vector<nolife_planet *>;
  vector<nointelligence_planet *> *Nointelligence_planets = new vector<nointelligence_planet *>;
  vector<intelligence_planet *>   *Intelligence_planets   = new vector<intelligence_planet *>;
  vector<SNe *>                   *SNe_list               = new vector<SNe *>;
  vector<SNe *>                   *excute_list            = new vector<SNe *>;
  // vector<vector<SNe *> *> *        total_SNe_list         = new vector<vector<SNe *> *>;
};

struct location_list {
  location *p_location_list = new location[160000];
  int       location_count  = 0;
};

void galaxy_formation(vector<location> *Location_List);