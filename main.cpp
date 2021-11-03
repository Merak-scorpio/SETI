#include "SNe.hpp"
#include "constance.hpp"
#include "galaxy.hpp"
#include "planet.hpp"


#include <omp.h>

using namespace std;


default_random_engine             e_main(time(0));
uniform_real_distribution<double> u_SNe_num(0, SNe_list_count - 1);

int main() {

  constance_set();

  cout << "T_evo= " << T_evo << endl;
  cout << "lambda_A= " << lambda_A << endl;
  cout << "P_ann= " << P_ann << endl;

  double   start = omp_get_wtime();
  ofstream ofs;
  string   str_year;

  // 1   初始化银河系
  cout << "initialize galaxy" << endl;
  vector<location> Galaxy;
  galaxy_formation(&Galaxy);
  double t_galaxy_formation = omp_get_wtime();
  cout << "initialization cost " << t_galaxy_formation - start << " s" << endl;


  cout << "Generate supernova initial data" << endl;
  vector<vector<SNe *> *> total_SNe_list;
  vector<int>             SNe_num(evo_time);
  for (int year = 1; year <= evo_time; year++) {
    SNe_num[year - 1] = u_SNe_num(e_main);
  }
  for (int i = 0; i < SNe_list_count; i++) {
    vector<SNe *> *SNe_list = new vector<SNe *>;
    SNe_loop(SNe_list);
    total_SNe_list.push_back(SNe_list);
  }
  double t_SNe_formation = omp_get_wtime();
  cout << "Supernova initial data generation completed, cost " << t_SNe_formation - t_galaxy_formation << " s" << endl;


  cout << "Record supernova data" << endl;
#pragma omp parallel for
  for (int i = 0; i < int(Galaxy.size()); i++) {
    for (int count = 0; count < SNe_list_count; count++) {
      SNe_add(&Galaxy[i], total_SNe_list[count]);
    }
  }
  double t_SNe_record = omp_get_wtime();
  cout << "Supernova recording completed, cost " << t_SNe_record - t_SNe_formation << " s" << endl;


  // 2   演化开始(每一步时长 1Myr)
  for (int year = 1; year <= evo_time; year++) {
    cout << "No. " << year << " million years start" << endl;

    if (year <= 10000) {
      // 2.1 产生新的恒星和行星
      cout << "Producing planets" << endl;
      double t_planet_start = omp_get_wtime();
#pragma omp parallel for
      for (int i = 0; i < int(Galaxy.size()); i++) {
        planet_loop(&Galaxy[i]);
      }
      double t_planet_finish = omp_get_wtime();
      cout << "Producted, cost " << t_planet_finish - t_planet_start << " s" << endl;
    }

    // 2.2 生命进程判定
    double t_process_start = omp_get_wtime();
    cout << "Life processing" << endl;
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      planet_nolife_process(&Galaxy[i]);
      planet_nointelligence_process(&Galaxy[i]);
      planet_intelligence_process(&Galaxy[i]);
    }
    double t_process_finish = omp_get_wtime();
    cout << "Process finished, cost " << t_process_finish - t_process_start << " s" << endl;


    // 2.4 超新星爆炸，重置灭绝半径内行星上的生命进程
    double t_SNe_start = omp_get_wtime();
    cout << "SNe excuting" << endl;
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      SNe_excute(&Galaxy[i], SNe_num[year - 1]);
    }
    double t_SNe_finish = omp_get_wtime();
    cout << "Excuted, cost" << t_SNe_finish - t_SNe_start << " s" << endl;


    // 2.5 记录演化结果
    if (year % 1 == 0) {
      str_year = to_string(year);
      ofs.open("./result/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A + "/P_ann__" + str_P_ann + "/result_" + str_year + ".csv");
      ofs << "x,y,age" << endl;
      for (int i = 0; i < int(Galaxy.size()); i++) {
        if (!Galaxy[i].Intelligence_planets->empty()) {
          for (int j = 0; j < int(Galaxy[i].Intelligence_planets->size()); j++) {
            ofs << (*Galaxy[i].Intelligence_planets)[j]->x_coordinate << ',' << (*Galaxy[i].Intelligence_planets)[j]->y_coordinate << ','
                << (*Galaxy[i].Intelligence_planets)[j]->intelligence_age << endl;
          }
        }
      }
      ofs.close();
    }
  }


  double finish = omp_get_wtime();
  cout << "Total cost" << finish - start << " s" << endl;

  return 0;
}