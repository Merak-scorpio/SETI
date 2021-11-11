#include "SNe.hpp"

default_random_engine e_SNe(time(0));

uniform_real_distribution<double> u_(0, sigma_c);
uniform_real_distribution<double> u_r(0, R);
uniform_real_distribution<double> u_theta(-M_PI, M_PI);
uniform_real_distribution<double> u_z(0, thick);
uniform_real_distribution<double> u_10_100_(10, 100);
normal_distribution<double>       n(M_SNII, 1.35);

void SNe_loop(vector<SNe *> *SNe_list) {

  double r, theta, temp_judge, f_r_theta;

  int count = 0;
  while (count < SNe_total) {

    r          = u_r(e_SNe);
    theta      = u_theta(e_SNe);
    temp_judge = u_(e_SNe);
    f_r_theta  = sigma_c * exp(-r / h_R);

    // if (temp_judge < f_r_theta) {
    //   SNe *temp_SNe          = new SNe;
    //   temp_SNe->x_coordinate = r * cos(theta);
    //   temp_SNe->y_coordinate = r * sin(theta);
    //   temp_SNe->z_coordinate = u_z(e_SNe);
    //   temp_SNe->D_SNe        = d_SNII * exp(-0.4 * (n(e_SNe) - M_SNII));
    //   SNe_list->emplace_back(temp_SNe);
    //   count++;
    // }

    if (temp_judge < f_r_theta) {
      (*SNe_list)[count]->x_coordinate = r * cos(theta);
      (*SNe_list)[count]->y_coordinate = r * sin(theta);
      (*SNe_list)[count]->z_coordinate = u_z(e_SNe);
      (*SNe_list)[count]->D_SNe        = d_SNII * exp(-0.4 * (n(e_SNe) - M_SNII));
      count++;
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
  // vector<SNe *> *temp_SNe_list = new vector<SNe *>;
  for (int i = 0; i < SNe_total; i++) {
    if (SNe_add_judge(location, (*SNe_list)[i])) {
      SNe *temp_SNe          = new SNe;
      temp_SNe->x_coordinate = (*SNe_list)[i]->x_coordinate;
      temp_SNe->y_coordinate = (*SNe_list)[i]->y_coordinate;
      temp_SNe->z_coordinate = (*SNe_list)[i]->z_coordinate;
      temp_SNe->D_SNe        = (*SNe_list)[i]->D_SNe;
      location->SNe_list->emplace_back(temp_SNe);
    }
  }
  // location->total_SNe_list->emplace_back(temp_SNe_list);
}

bool excute_judge(double x1, double y1, double z1, double x2, double y2, double z2, double r) {
  double d_2 = (x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2);
  return d_2 < r * r;
}

// void SNe_excute(location *location, int SNe_num) {
//   vector<int> Nointelligence_erase_list;
//   vector<int> Intelligence_erase_list;

//   for (int i = 0; i < int(location->Nolife_planets->size()); i++) {
//     for (int j = 0; j < int((*location->total_SNe_list)[SNe_num]->size()); j++) {
//       if (excute_judge((*location->Nolife_planets)[i]->x_coordinate, (*location->Nolife_planets)[i]->y_coordinate,
//                        (*location->Nolife_planets)[i]->z_coordinate, (*(*location->total_SNe_list)[SNe_num])[j]->x_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->y_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->z_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->D_SNe)) {

//         (*location->Nolife_planets)[i]->planet_T_min  = u_10_100_(e_SNe);
//         (*location->Nolife_planets)[i]->life_evo_time = 0;
//         break;
//       }
//     }
//   }

//   for (int i = 0; i < int(location->Nointelligence_planets->size()); i++) {
//     for (int j = 0; j < int((*location->total_SNe_list)[SNe_num]->size()); j++) {
//       if (excute_judge((*location->Nointelligence_planets)[i]->x_coordinate, (*location->Nointelligence_planets)[i]->y_coordinate,
//                        (*location->Nointelligence_planets)[i]->z_coordinate, (*(*location->total_SNe_list)[SNe_num])[j]->x_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->y_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->z_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->D_SNe)) {

//         Nointelligence_erase_list.emplace_back(i);
//         nolife_planet *temp_planet = new nolife_planet;
//         temp_planet->planet_T_min  = u_10_100_(e_SNe);
//         temp_planet->star_age      = (*location->Nointelligence_planets)[i]->star_age;
//         temp_planet->star_T_L      = (*location->Nointelligence_planets)[i]->star_T_L;
//         temp_planet->x_coordinate  = (*location->Nointelligence_planets)[i]->x_coordinate;
//         temp_planet->y_coordinate  = (*location->Nointelligence_planets)[i]->y_coordinate;
//         temp_planet->z_coordinate  = (*location->Nointelligence_planets)[i]->z_coordinate;
//         location->Nolife_planets->emplace_back(temp_planet);
//         break;
//       }
//     }
//   }

//   for (int i = 0; i < int(location->Intelligence_planets->size()); i++) {
//     for (int j = 0; j < int((*location->total_SNe_list)[SNe_num]->size()); j++) {
//       if (excute_judge((*location->Intelligence_planets)[i]->x_coordinate, (*location->Intelligence_planets)[i]->y_coordinate,
//                        (*location->Intelligence_planets)[i]->z_coordinate, (*(*location->total_SNe_list)[SNe_num])[j]->x_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->y_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->z_coordinate,
//                        (*(*location->total_SNe_list)[SNe_num])[j]->D_SNe)) {

//         Intelligence_erase_list.emplace_back(i);
//         nolife_planet *temp_planet = new nolife_planet;
//         temp_planet->planet_T_min  = u_10_100_(e_SNe);
//         temp_planet->star_age      = (*location->Intelligence_planets)[i]->star_age;
//         temp_planet->star_T_L      = (*location->Intelligence_planets)[i]->star_T_L;
//         temp_planet->x_coordinate  = (*location->Intelligence_planets)[i]->x_coordinate;
//         temp_planet->y_coordinate  = (*location->Intelligence_planets)[i]->y_coordinate;
//         temp_planet->z_coordinate  = (*location->Intelligence_planets)[i]->z_coordinate;
//         location->Nolife_planets->emplace_back(temp_planet);
//         break;
//       }
//     }
//   }

//   for (int i = 0; i < int(Nointelligence_erase_list.size()); i++) {

//     if ((*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i] != NULL) {
//       delete (*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i];
//       (*location->Nointelligence_planets)[Nointelligence_erase_list[i] - i] = NULL;
//     }
//     location->Nointelligence_planets->erase(location->Nointelligence_planets->begin() + Nointelligence_erase_list[i] - i);
//   }

//   for (int i = 0; i < int(Intelligence_erase_list.size()); i++) {
//     if ((*location->Intelligence_planets)[Intelligence_erase_list[i] - i] != NULL) {
//       delete (*location->Intelligence_planets)[Intelligence_erase_list[i] - i];
//       (*location->Intelligence_planets)[Intelligence_erase_list[i] - i] = NULL;
//     }
//     location->Intelligence_planets->erase(location->Intelligence_planets->begin() + Intelligence_erase_list[i] - i);
//   }
// }

void SNe_excute(location *location) {
  vector<int> Nointelligence_erase_list;
  vector<int> Intelligence_erase_list;

  for (int i = 0; i < int(location->Nolife_planets->size()); i++) {
    for (int j = 0; j < int(location->SNe_list->size()); j++) {
      if (excute_judge((*location->Nolife_planets)[i]->x_coordinate, (*location->Nolife_planets)[i]->y_coordinate,
                       (*location->Nolife_planets)[i]->z_coordinate, (*location->SNe_list)[j]->x_coordinate,
                       (*location->SNe_list)[j]->y_coordinate, (*location->SNe_list)[j]->z_coordinate, (*location->SNe_list)[j]->D_SNe)) {
        (*location->Nolife_planets)[i]->planet_T_min  = u_10_100_(e_SNe);
        (*location->Nolife_planets)[i]->life_evo_time = 0;
        break;
      }
    }
  }

  for (int i = 0; i < int(location->Nointelligence_planets->size()); i++) {
    for (int j = 0; j < int(location->SNe_list->size()); j++) {
      if (excute_judge((*location->Nointelligence_planets)[i]->x_coordinate, (*location->Nointelligence_planets)[i]->y_coordinate,
                       (*location->Nointelligence_planets)[i]->z_coordinate, (*location->SNe_list)[j]->x_coordinate,
                       (*location->SNe_list)[j]->y_coordinate, (*location->SNe_list)[j]->z_coordinate, (*location->SNe_list)[j]->D_SNe)) {

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
    for (int j = 0; j < int(location->SNe_list->size()); j++) {
      if (excute_judge((*location->Intelligence_planets)[i]->x_coordinate, (*location->Intelligence_planets)[i]->y_coordinate,
                       (*location->Intelligence_planets)[i]->z_coordinate, (*location->SNe_list)[j]->x_coordinate,
                       (*location->SNe_list)[j]->y_coordinate, (*location->SNe_list)[j]->z_coordinate, (*location->SNe_list)[j]->D_SNe)) {

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

  for (int i = 0; i < int(location->SNe_list->size()); i++) {
    if ((*location->SNe_list)[0] != NULL) {
      delete (*location->SNe_list)[0];
      (*location->SNe_list)[0] = NULL;
      location->SNe_list->erase(location->SNe_list->begin());
    }
  }
}