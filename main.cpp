#include "SNe.hpp"
#include "constance.hpp"
#include "galaxy.hpp"
#include "planet.hpp"

#include <fstream>
#include <omp.h>

using namespace std;


default_random_engine             e_main(time(0));
uniform_real_distribution<double> u_SNe_num(0, SNe_list_count - 1);

int main() {
  double   start = omp_get_wtime();
  ofstream ofs;
  string   str_year;

  // 1   初始化银河系
  cout << "正在初始化" << endl;
  vector<location> Galaxy;
  galaxy_formation(&Galaxy);
  double t_galaxy_formation = omp_get_wtime();
  cout << "银河系初始化完成，用时" << t_galaxy_formation - start << "秒" << endl;


  cout << "生成超新星初始数据" << endl;
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
  cout << "超新星初始数据生成完成，用时" << t_SNe_formation - t_galaxy_formation << "秒" << endl;


  cout << "记录超新星数据" << endl;
#pragma omp parallel for
  for (int i = 0; i < int(Galaxy.size()); i++) {
    for (int count = 0; count < SNe_list_count; count++) {
      SNe_add(&Galaxy[i], total_SNe_list[count]);
    }
  }
  double t_SNe_record = omp_get_wtime();
  cout << "超新星记录完成，用时" << t_SNe_record - t_SNe_formation << "秒" << endl;


  // 2   演化开始(每一步时长 1Myr)
  for (int year = 1; year <= evo_time; year++) {
    cout << "第 " << year << " 百万年开始" << endl;


    // 2.1 产生新的恒星和行星
    cout << "开始产生行星" << endl;
    double t_planet_start = omp_get_wtime();
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      planet_loop(&Galaxy[i]);
    }
    double t_planet_finish = omp_get_wtime();
    cout << "行星产生完毕，用时" << t_planet_finish - t_planet_start << "秒" << endl;


    // 2.2 生命进程判定
    double t_process_start = omp_get_wtime();
    cout << "开始判定行星生命进程" << endl;
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      planet_nolife_process(&Galaxy[i]);
      planet_nointelligence_process(&Galaxy[i]);
      planet_intelligence_process(&Galaxy[i]);
    }
    double t_process_finish = omp_get_wtime();
    cout << "生命进程判定完毕，用时" << t_process_finish - t_process_start << "秒" << endl;


    // 2.4 超新星爆炸，重置灭绝半径内行星上的生命进程
    double t_SNe_start = omp_get_wtime();
    cout << "开始超新星灭绝生命" << endl;
#pragma omp parallel for
    for (int i = 0; i < int(Galaxy.size()); i++) {
      SNe_excute(&Galaxy[i], SNe_num[year - 1]);
    }
    double t_SNe_finish = omp_get_wtime();
    cout << "超新星灭绝完毕，用时" << t_SNe_finish - t_SNe_start << "秒" << endl;


    // 2.5 记录演化结果
    if (year % 100 == 0) {
      str_year = to_string(year);
      ofs.open("./result/result_" + str_year + ".csv");
      ofs << "x,y" << endl;
      for (int i = 0; i < int(Galaxy.size()); i++) {
        if (!Galaxy[i].Intelligence_planets->empty()) {
          for (int j = 0; j < int(Galaxy[i].Intelligence_planets->size()); j++) {
            ofs << (*Galaxy[i].Intelligence_planets)[j]->x_coordinate << ',' << (*Galaxy[i].Intelligence_planets)[j]->y_coordinate << endl;
          }
        }
      }
      ofs.close();
    }
  }


  double finish = omp_get_wtime();
  cout << "总计用时" << finish - start << endl;
  return 0;
}
