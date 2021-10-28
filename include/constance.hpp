#pragma once

#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace std;

extern double cell_r;    //每个细分空间边长
extern double R;         //银河系半径(kpc)
extern double thick;     //银河系厚度
extern double V_cell;    //每个区域面积
extern int    A;         //恒星生成参数
extern double N;         //恒星生成参数
extern double sigma_c;   //
extern double h_R;       //
extern int    T_Lsun;    //太阳主序时间
extern double lambda_A;  //概率参数
extern int    T_evo;     //生命进化所需时间
extern double P_ann;     //每步自毁概率
extern int    evo_time;  //演化时长
extern int    T_step;    //演化步长
extern double P_Life;    //生命从无到有的泊松概率
extern int    SNe_total; //超新星总数
extern double d_SNII;    //超新星平均灭绝半径
extern double M_SNII;    //超新星评价绝对星等

extern int SNe_list_count; //超新星表容量大小