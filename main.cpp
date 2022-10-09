#include "SNe.hpp"
#include "constance.hpp"
#include "galaxy.hpp"
#include "planet.hpp"


#include <omp.h>

using namespace std;


default_random_engine e_main(time(0));

int main() {

  constance_set();
  cout << "T_evo= " << T_evo << endl;
  cout << "lambda_A= " << lambda_A << endl;
  cout << "P_ann= " << P_ann << endl;
  cout << "mass:" << str_star_mass_min << "--" << str_star_mass_max << endl;

  double   start = omp_get_wtime();
  ofstream total_ofs;
  ofstream ofs;
  total_ofs.open("./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A +
                     "/P_ann__" + str_P_ann + "/result.csv",
                 ios::app);
  total_ofs << "distence,year,count" << endl;
  total_ofs.close();
  string str_year;


  // 1   初始化银河系
  cout << "initialize galaxy" << endl;
  vector<location> Galaxy;
  vector<SNe *>    SNe_list;
  for (int i = 0; i < 25000; i++) {
    SNe *temp_SNe          = new SNe;
    temp_SNe->D_SNe        = 0;
    temp_SNe->x_coordinate = 0;
    temp_SNe->y_coordinate = 0;
    temp_SNe->z_coordinate = 0;
    SNe_list.emplace_back(temp_SNe);
  }
  galaxy_formation(&Galaxy);
  double t_galaxy_formation = omp_get_wtime();
  cout << "initialization cost " << t_galaxy_formation - start << " s" << endl;


  // 2   演化开始(每一步时长 1Myr)
  for (int year = 1; year <= evo_time; year++) {
    cout << "No. " << year << " million years start" << endl;

// 2.1 产生新的恒星和行星
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      planet_loop(&Galaxy[i], year);
    }

    // 2.2 生命进程判定
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      planet_nolife_process(&Galaxy[i]);
      planet_nointelligence_process(&Galaxy[i]);
      planet_intelligence_process(&Galaxy[i]);
    }


    str_year = to_string(year);

    /*
        if (year % 100 == 0) {
          ofstream planet_ofs;
          planet_ofs.open("./result/temp/planet/year_" + str_year + "_planet.csv");
          planet_ofs << "theta,d" << endl;
          for (int i = 0; i < int(Galaxy.size()); i++) {
            for (int j = 0; j < int(Galaxy[i].Nointelligence_planets->size()); j++) {
              double x     = (*Galaxy[i].Nointelligence_planets)[j]->x_coordinate;
              double y     = (*Galaxy[i].Nointelligence_planets)[j]->y_coordinate;
              double d     = sqrt(x * x + y * y);
              double theta = atan(x / y);
              if (x < 0) {
                theta = theta + M_PI;
              }
              theta = theta * 180 / M_PI;
              planet_ofs << theta << "," << d << endl;
            }
            for (int j = 0; j < int(Galaxy[i].Nolife_planets->size()); j++) {
              double x     = (*Galaxy[i].Nolife_planets)[j]->x_coordinate;
              double y     = (*Galaxy[i].Nolife_planets)[j]->y_coordinate;
              double d     = sqrt(x * x + y * y);
              double theta = atan(x / y);
              if (x < 0) {
                theta = theta + M_PI;
              }
              theta = theta * 180 / M_PI;
              planet_ofs << theta << "," << d << endl;
            }
            for (int j = 0; j < int(Galaxy[i].Intelligence_planets->size()); j++) {
              double x     = (*Galaxy[i].Intelligence_planets)[j]->x_coordinate;
              double y     = (*Galaxy[i].Intelligence_planets)[j]->y_coordinate;
              double d     = sqrt(x * x + y * y);
              double theta = atan(x / y);
              if (x < 0) {
                theta = theta + M_PI;
              }
              theta = theta * 180 / M_PI;
              planet_ofs << theta << "," << d << endl;
            }
          }
          planet_ofs.close();
        }
    */


    // 2.4 超新星爆炸，重置灭绝半径内行星上的生命进程

    SNe_loop(&SNe_list);
    /*
        if (year % 100 == 0) {
          ofstream SNe_ofs;
          SNe_ofs.open("./result/temp/SNe/year_" + str_year + "_SNe.csv");
          SNe_ofs << "theta,d" << endl;
          for (int i = 0; i < int(SNe_list.size()); i++) {

            double x     = SNe_list[i]->x_coordinate;
            double y     = SNe_list[i]->y_coordinate;
            double d     = sqrt(x * x + y * y);
            double theta = atan(x / y);
            if (x < 0) {
              theta = theta + M_PI;
            }
            theta = theta * 180 / M_PI;
            SNe_ofs << theta << "," << d << endl;
          }
          SNe_ofs.close();
        }
    */

#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      SNe_add(&Galaxy[i], &SNe_list);
      // cout << Galaxy[i].r << endl;
    }


#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      SNe_excute(&Galaxy[i]);
    }
    // SNe_end(&SNe_list, &Galaxy);

    // 2.5 记录演化结果
    if (year % 100 == 0) {
      vector<double> xaxis(bins);
      for (int i = 0; i < bins; i++) {
        xaxis[i] = i * 22.5 / bins;
      }
      vector<double> xcount(bins);
      for (int i = 0; i < bins; i++) {
        xcount[i] = 0;
      }
      int    local;
      double x, y, d, theta;
      str_year = to_string(year);
      ofs.open("./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A +
               "/P_ann__" + str_P_ann + "/result_" + str_year + ".csv");
      ofs << "theta,d,age" << endl;

      for (int i = 0; i < int(Galaxy.size()); i++) {
        if (!Galaxy[i].Intelligence_planets->empty()) {
          for (int j = 0; j < int(Galaxy[i].Intelligence_planets->size()); j++) {
            x     = (*Galaxy[i].Intelligence_planets)[j]->x_coordinate;
            y     = (*Galaxy[i].Intelligence_planets)[j]->y_coordinate;
            d     = sqrt(x * x + y * y);
            theta = atan(x / y);
            if (x < 0) {
              theta = theta + M_PI;
            }
            theta = theta * 180 / M_PI;
            ofs << theta << "," << d << "," << (*Galaxy[i].Intelligence_planets)[j]->intelligence_age << endl;
            local = int(d * bins / 22.5);
            xcount[local]++;
          }
        }
      }
      ofs.close();
      total_ofs.open("./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" +
                         str_lambda_A + "/P_ann__" + str_P_ann + "/_result_T_evo__" + str_T_evo + +"_lambda_A__" + str_lambda_A +
                         "_P_ann__" + str_P_ann + ".csv",
                     ios::app);
      for (int i = 0; i < bins; i++) {
        total_ofs << xaxis[i] << ',' << year << ',' << xcount[i] << endl;
      }
      total_ofs.close();
    }
  }

  double finish = omp_get_wtime();
  cout << "Total cost" << finish - start << " s" << endl;

  return 0;
}