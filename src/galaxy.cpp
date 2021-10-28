#include "galaxy.hpp"


void galaxy_formation(vector<location> *Location_List) {

  string      str_order;
  const char *order;


#ifdef WIN32
  str_order = "rd/s/q ./result/T_evo=" + str_T_evo + "/lambda_A=" + str_lambda_A + "/P_ann=" + str_P_ann;
  order     = str_order.data();
  system(order);
#elif linux
  str_order = "rm -rf ./result/T_evo=" + str_T_evo + "/lambda_A=" + str_lambda_A + "/P_ann=" + str_P_ann;
  order     = str_order.data();
  system(order);
#endif
  str_order = "mkdir ./result/T_evo=" + str_T_evo + "/lambda_A=" + str_lambda_A + "/P_ann=" + str_P_ann;
  order     = str_order.data();
  system(order);
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
