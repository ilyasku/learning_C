#include "harmonic_oscillator.h"

double V(double x){
  return -0.5 * K * x * x;
}

double psi_n_plus(double xi_n_minus, double psi_n_minus,
		  double xi_n, double psi_n,
		  double xi_n_plus){
  double g_n_plus = -2 * epsilon + xi_n_plus * xi_n_plus;
  double f_n_plus = 1 + g_n_plus * delta_xi_squared / 12.;
  double g_n = -2 * epsilon + xi_n * xi_n;
  double f_n = 1 + g_n * delta_xi_squared / 12.;
  double g_n_minus = -2 * epsilon + xi_n_minus * xi_n_minus;
  double f_n_minus = 1 + g_n_minus * delta_xi_squared / 12.;
  return ((12 - 10 * f_n) * psi_n - f_n_minus * psi_n_minus) / f_n_plus;
}
