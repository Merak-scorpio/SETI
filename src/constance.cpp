#include "constance.hpp"

double cell_r         = 0.1;                  //每个细分空间边长
double R              = 22.5;                 //银河系半径(kpc)
double thick          = 0.7;                  //银河系厚度
double V_cell         = cell_r * cell_r;      //每个区域面积
int    A              = 250;                  //恒星生成参数
double N              = 1.4;                  //恒星生成参数
double sigma_c        = 1.4 * pow(10, 2);     //
double h_R            = 2.25;                 //
int    T_Lsun         = 11000;                //太阳主序时间
double lambda_A       = 1;                    //概率参数
int    T_evo          = 20;                   //生命进化所需时间
double P_ann          = 0.0001;               //每步自毁概率
int    evo_time       = 100;                  //演化时长
int    T_step         = 1;                    //演化步长
double P_Life         = 1 - exp(-lambda_A);   //生命从无到有的泊松概率
double d_SNII         = 0.08;                 //超新星平均灭绝半径
double M_SNII         = -16.89;               //超新星评价绝对星等
int    SNe_total      = int(0.025 * 1000000); //超新星总数
int    SNe_list_count = 50;                   //超新星表容量大小