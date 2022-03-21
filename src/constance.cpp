#include "constance.hpp"

double cell_r    = 0.1;                  //每个细分空间边长
double R         = 22.5;                 //银河系半径(kpc)
double core_R    = 2.5;                  //核球半径
double thick     = 0.7;                  //银河系厚度
double V_cell    = cell_r * cell_r;      //每个区域面积
double A         = 2.5 * pow(10, -6.4);  //恒星生成参数
double N         = 1.4;                  //恒星生成参数
double sigma_c   = 1.4 * pow(10, 2);     //恒星分布参数
double h_R       = 2.25;                 //银河系标长
int    T_Lsun    = 11000;                //太阳主序时间
double d_SNII    = 0.08;                 //超新星平均灭绝半径
double M_SNII    = -16.89;               //超新星评价绝对星等
int    SNe_total = int(0.025 * 1000000); //超新星总数
double k         = 2.5 * 1000;           //内落时标

int    T_evo;          //生命进化所需时间
double lambda_A;       //概率参数
double P_ann;          //每步自毁概率
int    evo_time;       //演化时长
double P_Life;         //生命从无到有的泊松概率
int    SNe_list_count; //超新星表容量大小
int    star_end_year;  //产生新恒星停止时间

int bins = 30; //画图精细度

string str_T_evo;
string str_lambda_A;
string str_P_ann;
string str_SNe_list_count;
string str_evo_time;
string str_star_end_year;

void constance_set() {
  ifstream fin("constance.csv");
  string   line;
  getline(fin, line);
  getline(fin, line);
  istringstream  sin(line);
  vector<string> fields;
  string         field;
  while (getline(sin, field, ',')) {
    fields.push_back(field);
  }

  str_T_evo    = fields[0];
  str_lambda_A = fields[1];
  str_P_ann    = fields[2];
  str_evo_time = fields[3];

  T_evo    = stoi(str_T_evo);
  lambda_A = stod(str_lambda_A);
  P_ann    = stod(str_P_ann);
  evo_time = stoi(str_evo_time);


  P_Life       = 1 - exp(-lambda_A);
  str_T_evo    = to_string(T_evo);
  str_lambda_A = to_string(lambda_A);
  str_P_ann    = to_string(P_ann);
}