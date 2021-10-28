#include "galaxy.hpp"


void galaxy_formation(vector<location> *Location_List) {

#ifdef WIN32
  system("rd/s/q .\\result");
#elif linux
  system("rm -rf result");
  system("mkdir result");
#endif

  cout << "清理完了" << endl;

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
