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
    fprintf(stderr, "last change of sign too far!\n");
    fprintf(stderr, "index = %i\n", index_of_last_sign_change);
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

int compute_psi(int N_intervals,		 
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

  return compute_psi(N_intervals, n_max_iterations - 1, expected_number_of_nodes);  
}

int find_classical_limit(int lower_limit, int upper_limit) {
  int index_found = 0;
  while (index_found == 0) {
    if ((upper_limit - lower_limit) <= 1) {
      return upper_limit;
    }    
    int index_classical_limit = (upper_limit + lower_limit) / 2;
    double v = vpot[index_classical_limit];
    if (v == epsilon) {
      return index_classical_limit;
    }
    if (v < epsilon) {
      lower_limit = index_classical_limit;
      continue;
    }
    upper_limit = index_classical_limit;    
  }
  fprintf(stderr,							\
	  "*** '%s' in '%s' on line %d failed with error '%s'.\n",	\
	  "ERROR", __FILE__, __LINE__,					\
	  "Could not find the index where V > E!");			\
  abort();
}

int compute_psi_from_right_to_left(int N_intervals, int index_classical_limit,
				   int n_max_iterations, double threshold_difference_in_derivatives) {

  int N_intervals_right = N_intervals - index_classical_limit;
  
  for (int i = 2; i >= N_intervals_right; i++) {
    psi_right[N_intervals_right - i] =
      ((12. - 10. * f[N_intervals - i + 1]) * psi_right[N_intervals_right - i + 1] - \
       psi_right[N_intervals_right - i + 2] * f[N_intervals + 2]) / f[N_intervals - i];
  }

  psi_matching_point_right = psi_right[0];
  double matching_factor = psi_matching_point_left / psi_matching_point_right;

  psi_right[0] = psi_matching_point_left;
  psi_right[1] = psi_right[1] * matching_factor;

  double difference_in_left_and_right_derivative =
    (psi[index_classical_limit - 1] + psi_right[1] \
     - (14. - 12. * f[index_classical_limit]) * psi_matching_point_left) / delta_xi;

  if (abs(difference_in_left_and_right_derivative) <= threshold_difference_in_derivatives ||
      n_max_iterations <= 1){
    normalize_psi_right_and_write_to_psi(matching_factor, N_intervals_right,
					 index_classical_limit);
    return n_max_iterations - 1;
  }

  if (difference_in_left_and_right_derivative < 0) {
    epsilon_min = epsilon;
  } else {
    epsilon_max = epsilon;
  }

  epsilon = 0.5 * (epsilon_max + epsilon_min);

  return compute_psi_from_right_to_left(N_intervals, index_classical_limit,
					n_max_iterations - 1, threshold_difference_in_derivatives);    
}

void normalize_psi_right_and_write_to_psi(double matching_factor, int N_intervals_right,
					  int index_classical_limit) {
  for (int i = 1; i <= N_intervals_right; i++) {
    psi[index_classical_limit + i] = psi_right[i] * matching_factor;
  }  
}

void normalize_psi(int N_intervals) {
  double norm = 0.;
  for (int i = 1; i <= N_intervals; i++) {
    norm += psi[i] * psi[i];
  }
  norm = delta_xi * (2. * norm + psi[0] * psi[0]);
  norm = sqrt(norm);

  for (int i = 0; i <= N_intervals; i++) {
    psi[i] /= norm;
  }
}
