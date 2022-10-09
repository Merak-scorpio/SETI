#include "galaxy.hpp"


void galaxy_formation(vector<location> *Location_List) {
  ofstream ofs;

  string      str_dirname;
  string      str_order;
  const char *dirname;
  const char *order;


#ifdef WIN32
  str_dirname = ".\\result";
  dirname     = str_dirname.data();
  if (_access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = ".\\result\\T_evo__" + str_T_evo;
  dirname     = str_dirname.data();
  if (_access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = ".\\result\\T_evo__" + str_T_evo + "\\lambda_A__" + str_lambda_A;
  dirname     = str_dirname.data();
  if (_access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = ".\\result\\T_evo__" + str_T_evo + "\\lambda_A__" + str_lambda_A + "\\P_ann__" + str_P_ann;
  dirname     = str_dirname.data();
  if (_access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_order = "rd/s/q " + str_dirname;
  order     = str_order.data();
  system(order);
  str_order = "mkdir .\\result\\T_evo__" + str_T_evo + "\\lambda_A__" + str_lambda_A + "\\P_ann__" + str_P_ann;
  order     = str_order.data();
  system(order);

  ofs.open(".\\result\\T_evo__" + str_T_evo + "\\lambda_A__" + str_lambda_A + "\\P_ann__" + str_P_ann + "\\_result_T_evo__" + str_T_evo +
           +"_lambda_A__" + str_lambda_A + "_P_ann__" + str_P_ann + ".csv");
  ofs << "distance,year,count" << endl;
  ofs.close();
#elif linux
  str_dirname = "./result";
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/mass_" + str_star_mass_min + "--" + str_star_mass_max;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A +
                "/P_ann__" + str_P_ann;
  dirname = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_order = "rm -rf " + str_dirname;
  order     = str_order.data();
  system(order);
  str_order = "mkdir ./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" +
              str_lambda_A + "/P_ann__" + str_P_ann;
  order = str_order.data();
  system(order);

  ofs.open("./result/mass_" + str_star_mass_min + "--" + str_star_mass_max + "/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A +
           "/P_ann__" + str_P_ann + "/_result_T_evo__" + str_T_evo + +"_lambda_A__" + str_lambda_A + "_P_ann__" + str_P_ann + ".csv");
  ofs << "distance,year,count" << endl;
  ofs.close();
#endif


  int    r2, n, x, y, size;
  double r;

  n    = int(R / cell_r);
  size = 2 * n + 1;
  for (int i = 0; i < size; i++) {
    x = i - n;
    for (int j = 0; j < size; j++) {
      r2 = (i - n) * (i - n) + (j - n) * (j - n);
      r  = sqrt(r2) * cell_r;
      if (r < R && r > 0) {
        y = j - n;

        location temp_location;
        temp_location.x   = x;
        temp_location.y   = y;
        temp_location.r   = r;
        temp_location.SFR = 0;
        temp_location.Normal_A =
            (3.7 * pow(10, 7) * exp((8 - r) / h_R) + sigma_c * exp(-r / h_R)) / ((k * r / h_R) * (1 - exp(-13500 * h_R / (k * r))));

        Location_List->push_back(temp_location);
      }
    }
  }
}
