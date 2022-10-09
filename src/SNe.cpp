#include "SNe.hpp"

default_random_engine             e_SNe(time(0));
uniform_real_distribution<double> u_(0, sigma_c);
uniform_real_distribution<double> u_r(0, R);
uniform_real_distribution<double> u_theta(-M_PI, M_PI);
uniform_real_distribution<double> u_z(0, thick);
uniform_real_distribution<double> u_10_100_(10, 100);
normal_distribution<double>       n(M_SNII, 1.35);


void SNe_loop(vector<SNe *> *SNe_list, vector<location> *Galaxy) {
  for (int i = 0; i < int(Galaxy->size()); i++) {
    vector<int> large_star_erase_list;
    for (int j = 0; j < int((*Galaxy)[i].Large_mass_stars->size()); j++) {
      (*(*Galaxy)[i].Large_mass_stars)[j]->star_age++;
      if ((*(*Galaxy)[i].Large_mass_stars)[j]->star_age > (*(*Galaxy)[i].Large_mass_stars)[j]->star_T_L) {
        SNe *temp_SNe          = new SNe;
        temp_SNe->x_coordinate = (*(*Galaxy)[i].Large_mass_stars)[j]->x_coordinate;
        temp_SNe->y_coordinate = (*(*Galaxy)[i].Large_mass_stars)[j]->y_coordinate;
        temp_SNe->z_coordinate = (*(*Galaxy)[i].Large_mass_stars)[j]->z_coordinate;
        temp_SNe->D_SNe        = d_SNII * pow(pow(10, -0.4 * (n(e_SNe) - M_SNII)), 0.5);
        SNe_list->emplace_back(temp_SNe);
        large_star_erase_list.emplace_back(j);
      }
    }

    for (int j = 0; j < int(large_star_erase_list.size()); j++) {
      if ((*(*Galaxy)[i].Large_mass_stars)[large_star_erase_list[j] - j] != NULL) {
        delete (*(*Galaxy)[i].Large_mass_stars)[large_star_erase_list[j] - j];
        (*(*Galaxy)[i].Large_mass_stars)[large_star_erase_list[j] - j] = NULL;
      }
      (*Galaxy)[i].Large_mass_stars->erase((*Galaxy)[i].Large_mass_stars->begin() + large_star_erase_list[j] - j);
    }
  }
}

bool SNe_add_judge(location *location, SNe *SNe) {
  double cir_x_coordinate = SNe->x_coordinate;
  double cir_y_coordinate = SNe->y_coordinate;
  double cir_r            = SNe->D_SNe;

  double sqr_x_coordinate = location->x * cell_r;
  double sqr_y_coordinate = location->y * cell_r;
  double half_sqr_l       = 0.5 * cell_r;

  double v[2]   = {abs(cir_x_coordinate - sqr_x_coordinate), abs(cir_y_coordinate - sqr_y_coordinate)};
  double u[2]   = {max<double>((v[0] - half_sqr_l), 0), max<double>((v[1] - half_sqr_l), 0)};
  double u_len2 = (u[0] * u[0]) + (u[1] * u[1]);

  return u_len2 <= (cir_r * cir_r);
}

void SNe_add(location *location, vector<SNe *> *SNe_list) {
  for (int i = 0; i < int(SNe_list->size()); i++) {
    if (SNe_add_judge(location, (*SNe_list)[i])) {
      SNe *temp_SNe          = new SNe;
      temp_SNe->x_coordinate = (*SNe_list)[i]->x_coordinate;
      temp_SNe->y_coordinate = (*SNe_list)[i]->y_coordinate;
      temp_SNe->z_coordinate = (*SNe_list)[i]->z_coordinate;
      temp_SNe->D_SNe        = (*SNe_list)[i]->D_SNe;
      location->excute_list->emplace_back(temp_SNe);
    }
  }
}

bool excute_judge(double x1, double y1, double z1, double x2, double y2, double z2, double r) {
  double d_2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
  return d_2 < r * r;
}

void SNe_excute(location *location) {

  vector<int> Nointelligence_erase_list;
  vector<int> Intelligence_erase_list;

  for (int i = 0; i < int(location->Nolife_planets->size()); i++) {
    for (int j = 0; j < int(location->excute_list->size()); j++) {
      if (excute_judge((*location->Nolife_planets)[i]->x_coordinate, (*location->Nolife_planets)[i]->y_coordinate,
                       (*location->Nolife_planets)[i]->z_coordinate, (*location->excute_list)[j]->x_coordinate,
                       (*location->excute_list)[j]->y_coordinate, (*location->excute_list)[j]->z_coordinate,
                       (*location->excute_list)[j]->D_SNe)) {
        (*location->Nolife_planets)[i]->planet_T_min  = u_10_100_(e_SNe);
        (*location->Nolife_planets)[i]->life_evo_time = 0;
        break;
      }
    }
  }

  for (int i = 0; i < int(location->Nointelligence_planets->size()); i++) {
    for (int j = 0; j < int(location->excute_list->size()); j++) {
      if (excute_judge((*location->Nointelligence_planets)[i]->x_coordinate, (*location->Nointelligence_planets)[i]->y_coordinate,
                       (*location->Nointelligence_planets)[i]->z_coordinate, (*location->excute_list)[j]->x_coordinate,
                       (*location->excute_list)[j]->y_coordinate, (*location->excute_list)[j]->z_coordinate,
                       (*location->excute_list)[j]->D_SNe)) {
        Nointelligence_erase_list.emplace_back(i);
        nolife_planet *temp_planet = new nolife_planet;
        temp_planet->planet_T_min  = u_10_100_(e_SNe);
        temp_planet->star_age      = (*location->Nointelligence_planets)[i]->star_age;
        temp_planet->star_T_L      = (*location->Nointelligence_planets)[i]->star_T_L;
        temp_planet->x_coordinate  = (*location->Nointelligence_planets)[i]->x_coordinate;
        temp_planet->y_coordinate  = (*location->Nointelligence_planets)[i]->y_coordinate;
        temp_planet->z_coordinate  = (*location->Nointelligence_planets)[i]->z_coordinate;
        location->Nolife_planets->emplace_back(temp_planet);
        break;
      }
    }
  }

  for (int i = 0; i < int(location->Intelligence_planets->size()); i++) {
    for (int j = 0; j < int(location->excute_list->size()); j++) {
      if (excute_judge((*location->Intelligence_planets)[i]->x_coordinate, (*location->Intelligence_planets)[i]->y_coordinate,
                       (*location->Intelligence_planets)[i]->z_coordinate, (*location->excute_list)[j]->x_coordinate,
                       (*location->excute_list)[j]->y_coordinate, (*location->excute_list)[j]->z_coordinate,
                       (*location->excute_list)[j]->D_SNe)) {

        Intelligence_erase_list.emplace_back(i);
        nolife_planet *temp_planet = new nolife_planet;
        temp_planet->planet_T_min  = u_10_100_(e_SNe);
        temp_planet->star_age      = (*location->Intelligence_planets)[i]->star_age;
        temp_planet->star_T_L      = (*location->Intelligence_planets)[i]->star_T_L;
        temp_planet->x_coordinate  = (*location->Intelligence_planets)[i]->x_coordinate;
        temp_planet->y_coordinate  = (*location->Intelligence_planets)[i]->y_coordinate;
        temp_planet->z_coordinate  = (*location->Intelligence_planets)[i]->z_coordinate;
        location->Nolife_planets->emplace_back(temp_planet);
        break;
      }
    }
  }

  for (int i = 0; i < int(Nointelligence_erase_list.size()); i++) {

    if ((*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i] != NULL) {
      delete (*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i];
      (*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i] = NULL;
    }
    location->Nointelligence_planets->erase(location->Nointelligence_planets->begin() + Nointelligence_erase_list[i] - i);
  }

  for (int i = 0; i < int(Intelligence_erase_list.size()); i++) {
    if ((*location->Intelligence_planets)[Intelligence_erase_list[i] - i] != NULL) {
      delete (*location->Intelligence_planets)[Intelligence_erase_list[i] - i];
      (*location->Intelligence_planets)[Intelligence_erase_list[i] - i] = NULL;
    }
    location->Intelligence_planets->erase(location->Intelligence_planets->begin() + Intelligence_erase_list[i] - i);
  }

  for (int i = 0; i < int(location->excute_list->size()); i++) {
    if ((*location->excute_list)[0] != NULL) {
      delete (*location->excute_list)[0];
      (*location->excute_list)[0] = NULL;
      location->excute_list->erase(location->excute_list->begin());
    }
  }
}

void SNe_end(vector<SNe *> *SNe_list, vector<location> *Galaxy) {
  if (int(SNe_list->size()) != 0) {
    for (int i = 0; i < int(Galaxy->size()); i++) {
      for (int j = 0; j < int((*Galaxy)[i].SNe_list->size()); j++) {
        if ((*(*Galaxy)[i].SNe_list)[0] != NULL) {
          (*(*Galaxy)[i].SNe_list)[0] = NULL;
          (*Galaxy)[i].SNe_list->erase((*Galaxy)[i].SNe_list->begin());
        }
      }
    }
  }
}
