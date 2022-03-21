#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>
#include <stdio.h>
#include <string>
#include <vector>


using namespace std;

extern double cell_r;    //每个细分空间边长
extern double R;         //银河系半径(kpc)
extern double core_R;    //核球半径
extern double thick;     //银河系厚度
extern double h_R;       //银河系标长
extern double k;         //内落时标
extern double V_cell;    //每个区域面积
extern double A;         //恒星生成参数
extern double N;         //恒星生成参数
extern double sigma_c;   //恒星分布参数
extern int    T_Lsun;    //太阳主序时间
extern double P_Life;    //生命从无到有的泊松概率
extern int    SNe_total; //超新星总数
extern double d_SNII;    //超新星平均灭绝半径
extern double M_SNII;    //超新星评价绝对星等

extern double lambda_A; //概率参数
extern int    T_evo;    //生命进化所需时间
extern double P_ann;    //每步自毁概率
extern int    evo_time; //演化时长

extern int bins; //画图精细度

extern string str_T_evo;
extern string str_lambda_A;
extern string str_P_ann;

void constance_set();