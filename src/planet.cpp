#include "planet.hpp"


using namespace std;

default_random_engine e(time(0));

uniform_real_distribution<double> u_planet_phi(0.000448, 2.0939);
uniform_real_distribution<double> u_0_1(0, 1);
uniform_real_distribution<double> u_10_100(10, 100);
normal_distribution<double>       n_(M_SNII, 1.35);

double IMF_rev(double phi) {
  double m = 0;
  if (2.094 >= phi && phi >= 0.578) {
    m = pow((phi / 0.332), (-1 / 0.8));
  } else if (0.578 > phi && phi >= 0.178) {
    m = pow((phi / 0.178), (-1 / 1.7));
  } else if (0.178 > phi && phi >= 0.000447) {
    m = pow((phi / 0.178), (-1 / 1.3));
  }
  return m;
}


void planet_loop(location *location, int year) {


  double planet_x_coordinate, planet_y_coordinate, planet_z_coordinate, sne_x_coordinate, sne_y_coordinate, sne_z_coordinate, star_T_L,
      planet_T_min, temp_random, fin_t;

  fin_t = location->Normal_A * exp(-year * h_R / (k * location->r));
  location->Gas_mass_left += fin_t * V_cell;
  location->sigma_Gas = location->Gas_mass_left / V_cell;
  location->sigma_SFR = A * pow(location->sigma_Gas, 1.4) * V_cell + location->SFR_left;
  while (true) {
    double temp_star_mass_phi = u_planet_phi(e);
    double temp_star_mass     = IMF_rev(temp_star_mass_phi);
    if (location->sigma_SFR < temp_star_mass) {
      location->SFR_left = location->sigma_SFR;
      break;
    }
    location->sigma_SFR -= temp_star_mass;
    location->Gas_mass_left -= temp_star_mass;
    temp_random = u_0_1(e);
    if (temp_star_mass > 0.8 && temp_star_mass < 1.2 && temp_random <= 0.00616) {
      planet_x_coordinate        = (location->x - 0.5 + u_0_1(e)) * cell_r;
      planet_y_coordinate        = (location->y - 0.5 + u_0_1(e)) * cell_r;
      planet_z_coordinate        = u_0_1(e) * thick;
      star_T_L                   = T_Lsun * pow(1 / temp_star_mass, 2.5);
      planet_T_min               = u_10_100(e);
      nolife_planet *temp_planet = new nolife_planet;
      temp_planet->x_coordinate  = planet_x_coordinate;
      temp_planet->y_coordinate  = planet_y_coordinate;
      temp_planet->z_coordinate  = planet_z_coordinate;
      temp_planet->star_T_L      = star_T_L;
      temp_planet->planet_T_min  = planet_T_min;
      location->Nolife_planets->emplace_back(temp_planet);
    } else if (temp_star_mass > 8) {
      SNe *temp_sne          = new SNe;
      sne_x_coordinate       = (location->x - 0.5 + u_0_1(e)) * cell_r;
      sne_y_coordinate       = (location->y - 0.5 + u_0_1(e)) * cell_r;
      sne_z_coordinate       = u_0_1(e) * thick;
      temp_sne->x_coordinate = sne_x_coordinate;
      temp_sne->y_coordinate = sne_y_coordinate;
      temp_sne->z_coordinate = sne_z_coordinate;
      temp_sne->D_SNe        = d_SNII * exp(-0.4 * (n_(e) - M_SNII));
      location->SNe_list->emplace_back(temp_sne);
    }
  }
}

void planet_nolife_process(location *location) {
  double      temp_random;
  vector<int> erase_list;


  for (int i = 0; i < int(location->Nolife_planets->size()); i++) {
    (*location->Nolife_planets)[i]->star_age++;
    if ((*location->Nolife_planets)[i]->star_age > (*location->Nolife_planets)[i]->star_T_L) {
      erase_list.emplace_back(i);
      continue;
    }

    (*location->Nolife_planets)[i]->life_evo_time++;
    if ((*location->Nolife_planets)[i]->life_evo_time > (*location->Nolife_planets)[i]->planet_T_min) {
      temp_random = u_0_1(e);
      if (temp_random < P_Life) {

        nointelligence_planet *temp_planet = new nointelligence_planet;
        temp_planet->star_age              = (*location->Nolife_planets)[i]->star_age - 1;
        temp_planet->x_coordinate          = (*location->Nolife_planets)[i]->x_coordinate;
        temp_planet->y_coordinate          = (*location->Nolife_planets)[i]->y_coordinate;
        temp_planet->z_coordinate          = (*location->Nolife_planets)[i]->z_coordinate;
        temp_planet->star_T_L              = (*location->Nolife_planets)[i]->star_T_L;
        erase_list.emplace_back(i);
        location->Nointelligence_planets->emplace_back(temp_planet);
      }
    }
  }
  for (int i = 0; i < int(erase_list.size()); i++) {
    if ((*location->Nolife_planets)[erase_list[i] - i] != NULL) {
      delete (*location->Nolife_planets)[erase_list[i] - i];
      (*location->Nolife_planets)[erase_list[i] - i] = NULL;
    }
    location->Nolife_planets->erase(location->Nolife_planets->begin() + erase_list[i] - i);
  }
}

void planet_nointelligence_process(location *location) {
  vector<int> erase_list;

  for (int i = 0; i < int(location->Nointelligence_planets->size()); i++) {
    (*location->Nointelligence_planets)[i]->star_age++;
    if ((*location->Nointelligence_planets)[i]->star_age > (*location->Nointelligence_planets)[i]->star_T_L) {
      erase_list.emplace_back(i);
      continue;
    }

    (*location->Nointelligence_planets)[i]->life_age++;
    if ((*location->Nointelligence_planets)[i]->life_age > T_evo) {
      intelligence_planet *temp_planet = new intelligence_planet;
      temp_planet->star_age            = (*location->Nointelligence_planets)[i]->star_age - 1;
      temp_planet->x_coordinate        = (*location->Nointelligence_planets)[i]->x_coordinate;
      temp_planet->y_coordinate        = (*location->Nointelligence_planets)[i]->y_coordinate;
      temp_planet->z_coordinate        = (*location->Nointelligence_planets)[i]->z_coordinate;
      temp_planet->star_T_L            = (*location->Nointelligence_planets)[i]->star_T_L;
      erase_list.emplace_back(i);
      location->Intelligence_planets->emplace_back(temp_planet);
    }
  }
  for (int i = 0; i < int(erase_list.size()); i++) {
    if ((*location->Nointelligence_planets)[erase_list[i] - i] != NULL) {
      delete (*location->Nointelligence_planets)[erase_list[i] - i];
      (*location->Nointelligence_planets)[erase_list[i] - i] = NULL;
    }
    location->Nointelligence_planets->erase(location->Nointelligence_planets->begin() + erase_list[i] - i);
  }
}

void planet_intelligence_process(location *location) {
  double      temp_random;
  vector<int> erase_list;

  for (int i = 0; i < int(location->Intelligence_planets->size()); i++) {
    (*location->Intelligence_planets)[i]->star_age++;
    if ((*location->Intelligence_planets)[i]->star_age > (*location->Intelligence_planets)[i]->star_T_L) {
      erase_list.emplace_back(i);
      continue;
    }
    (*location->Intelligence_planets)[i]->intelligence_age++;
    temp_random = u_0_1(e);
    if (temp_random < P_ann) {
      erase_list.emplace_back(i);
    }
  }
  for (int i = 0; i < int(erase_list.size()); i++) {
    if ((*location->Intelligence_planets)[erase_list[i] - i] != NULL) {
      delete (*location->Intelligence_planets)[erase_list[i] - i];
      (*location->Intelligence_planets)[erase_list[i] - i] = NULL;
    }
    location->Intelligence_planets->erase(location->Intelligence_planets->begin() + erase_list[i] - i);
  }
}