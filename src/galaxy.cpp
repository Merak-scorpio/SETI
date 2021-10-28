#include "galaxy.hpp"


void galaxy_formation(vector<location> *Location_List) {
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
#elif linux
  str_dirname = "./result";
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/T_evo__" + str_T_evo;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_dirname = "./result/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A + "/P_ann__" + str_P_ann;
  dirname     = str_dirname.data();
  if (access(dirname, 0) != 0) {
    str_order = "mkdir " + str_dirname;
    order     = str_order.data();
    system(order);
  }
  str_order = "rm -rf " + str_dirname;
  order     = str_order.data();
  system(order);
  str_order = "mkdir ./result/T_evo__" + str_T_evo + "/lambda_A__" + str_lambda_A + "/P_ann__" + str_P_ann;
  order     = str_order.data();
  system(order);
#endif

  cout << "Initialization completed" << endl;

  int    r2, n, x, y, size;
  double r, sigma_Gas, sigma_SFR, Gas_Mass;

  n    = int(R / cell_r);
  size = 2 * n + 1;
  for (int i = 0; i < size; i++) {
    x = i - n;
    for (int j = 0; j < size; j++) {
      r2 = (i - n) * (i - n) + (j - n) * (j - n);
      r  = sqrt(r2) * cell_r;
      if (r < R) {
        y = j - n;

        sigma_Gas = sigma_c * exp(-r / h_R);
        sigma_SFR = A * pow(sigma_Gas, N);
        Gas_Mass  = sigma_SFR * V_cell;

        location temp_location;
        temp_location.x              = x;
        temp_location.y              = y;
        temp_location.each_step_mass = Gas_Mass;
        temp_location.Mass_left      = 0;

        Location_List->push_back(temp_location);
      }
    }
  }
}
