#include "harmonic_oscillator.h"

double get_g_n(double xi_n) {
  return -2 * epsilon + xi_n * xi_n;
}

double get_f_n(double g_n) {
  return 1 + g_n * delta_xi_squared / 12.;
}

double get_psi_1_for_nodes_even(double psi_0) {
  double g_0 = get_g_n(0.0);
  double g_1 = get_g_n(delta_xi);
  double f_0 = get_f_n(g_0);
  double f_1 = get_f_n(g_1);
  return (12 - 10 * f_0) * psi_0 / 2 / f_1;
}


double psi_n_plus(double xi_n_minus, double psi_n_minus,
		  double xi_n, double psi_n,
		  double xi_n_plus){
  double g_n_plus = get_g_n(xi_n_plus);
  double f_n_plus = get_f_n(g_n_plus);
  double g_n = get_g_n(xi_n);
  double f_n = get_f_n(g_n);
  double g_n_minus = get_g_n(xi_n_minus);
  double f_n_minus = get_f_n(g_n_minus);
  return ((12 - 10 * f_n) * psi_n - f_n_minus * psi_n_minus) / f_n_plus;
}

void compute_f(int N_intervals) {
  double ddx12 = delta_xi_squared / 12.;
  f[0] = ddx12 * (2. * (vpot[0] - epsilon));
  int index_of_last_sign_change = -1;
  for (int i = 1; i <= N_intervals; i++) {
    f[i] = ddx12 * 2. * (vpot[i] - epsilon);
    if (f[i] == 0) {
      f[i] = 1e-20;
    }
    if (f[i] != copysign(f[i], f[i - 1])) {
      index_of_last_sign_change = i;
    }
  }

  if (index_of_last_sign_change >= N_intervals - 2) {
    fprintf(stderr, "last change of sign too far!");
    exit(1);
  }

  if (index_of_last_sign_change < 1) {
    fprintf(stderr, "no classical turning point?");
    exit(1);
  }

  /* f(x) as required by the Numerov algorithm */
  for (int i = 0; i <= N_intervals; i++) {
    f[i] = 1. - f[i];
  }
}

int compute_psi(double epsilon_min, double epsilon_max,
		 int N_intervals,		 
		 int n_max_iterations,
		 int expected_number_of_nodes){

  epsilon = 0.5 * (epsilon_max + epsilon_min);

  compute_f(N_intervals);
  
  int count_sign_changes = 0;
  if (psi[0] == 0.0) {
    count_sign_changes += 1;
  } 
  
  for (int i = 1; i < N_intervals; i++) {
    psi[i + 1] = ((12. - f[i] * 10.) * psi[i] - f[i - 1] * psi[i - 1]) / f[i + 1];
    if (psi[i] != copysign(psi[i], psi[i + 1]))
      ++count_sign_changes;
  }

  if ((epsilon_max - epsilon_min) < 1e-10) {
    return n_max_iterations;
  }
  if (n_max_iterations == 1) {
    return n_max_iterations;
  }

  if (count_sign_changes > expected_number_of_nodes) {
    epsilon_max = epsilon;
  } else {
    epsilon_min = epsilon;
  }

  return compute_psi(epsilon_min, epsilon_max,
		     N_intervals, n_max_iterations - 1, expected_number_of_nodes);  
}
