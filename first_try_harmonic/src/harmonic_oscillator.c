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

void compute_f() {
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
    fprintf(stderr, "N_intervals = %i\n", N_intervals);
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

int compute_psi(int n_max_iterations,
		int expected_number_of_nodes,
		double threshold_in_derivatives)
{

  printf("--\n");
  printf("ENTER COMPUTE_PSI\n");
  printf("n_max_iterations = %i\n", n_max_iterations);

  
  epsilon = 0.5 * (epsilon_max + epsilon_min);

  find_classical_limit();

  if (index_classical_limit > N_intervals - 2) {
    fprintf(stderr, "ERROR: classical limit too close to xi_max\n");
    abort();
  }
  
  compute_f();
  
  int count_sign_changes = 0;
  if (psi[0] == 0.0) {
    count_sign_changes += 1;
  } 
  
  for (int i = 1; i < N_intervals; i++) {
    psi[i + 1] = ((12. - f[i] * 10.) * psi[i] - f[i - 1] * psi[i - 1]) / f[i + 1];
    if (psi[i] != copysign(psi[i], psi[i + 1]))
      ++count_sign_changes;
  }

  if (n_max_iterations <= 1) {
    compute_psi_right(DBL_MAX);
    concatenate_psi_left_and_psi_right();
    normalize_psi();    
    return 0;
  }
  
  if (count_sign_changes > expected_number_of_nodes) {
    epsilon_max = epsilon;
    return compute_psi(n_max_iterations - 1, expected_number_of_nodes,
		       threshold_in_derivatives);  
  }
  if (count_sign_changes < expected_number_of_nodes) {
    epsilon_min = epsilon;
    return compute_psi(n_max_iterations - 1, expected_number_of_nodes, threshold_in_derivatives);  
  }
      
  int psi_right_matches_psi_left = compute_psi_right(threshold_in_derivatives);

  if (psi_right_matches_psi_left == 0) {
    return compute_psi(n_max_iterations - 1, expected_number_of_nodes, threshold_in_derivatives);
  }
  
  concatenate_psi_left_and_psi_right();
  normalize_psi();
  return n_max_iterations;
}  

void concatenate_psi_left_and_psi_right(void) {
  for (int i = 1; i <= N_intervals - index_classical_limit; i++) {
    psi[index_classical_limit + i] = psi_right[i];
  }
}

void find_classical_limit() {
  int lower_limit = 0;
  int upper_limit = N_intervals;

  int index_found = 0;


  int count = 0;
  while (index_found == 0) {
    if ((upper_limit - lower_limit) <= 1) {
      fprintf(stderr, "return index %i after %i loops\n", upper_limit, count);
      index_classical_limit = upper_limit;
      return;
    }    
    index_classical_limit = (upper_limit + lower_limit) / 2;
    double v = vpot[index_classical_limit];
    if (v == epsilon) {
      fprintf(stderr, "due to v == epsilon,\nreturn index %i after %i loops\n", index_classical_limit, count);
      return;
    }
    if (v < epsilon) {
      lower_limit = index_classical_limit;
      continue;
    }
    upper_limit = index_classical_limit;
    count++;
  }
  fprintf(stderr,							\
	  "*** '%s' in '%s' on line %d failed with error '%s'.\n",	\
	  "ERROR", __FILE__, __LINE__,					\
	  "Could not find the index where V > E!");			\
  abort();
}

int compute_psi_right(double threshold_difference_in_derivatives) {

  printf("----\n");
  printf("ENTER COMPUTE_PSI_RIGHT\n");

  int N_intervals_right = N_intervals - index_classical_limit;

  printf("N_intervals_right = %i\n", N_intervals_right);
  
  if (PSI_RIGHT_ALLOCATED) free(psi_right);
  
  PSI_RIGHT_ALLOCATED = 1;
  psi_right = (double *) malloc((N_intervals_right + 1) * sizeof(double));
    
  psi_right[N_intervals_right] = delta_xi;
  psi_right[N_intervals_right - 1] =
    (12. - 10. * f[N_intervals]) * delta_xi / f[N_intervals - 1];
 
  psi_matching_point_left = psi[index_classical_limit];      
  
  for (int i = 2; i <= N_intervals_right; i++) {
    psi_right[N_intervals_right - i] =
      ((12. - 10. * f[N_intervals - i + 1]) * psi_right[N_intervals_right - i + 1] - \
       psi_right[N_intervals_right - i + 2] * f[N_intervals - i + 2]) / f[N_intervals - i];
  }
  
  psi_matching_point_right = psi_right[0];
  double matching_factor = psi_matching_point_left / psi_matching_point_right;

  psi_right[0] = psi_matching_point_left;
  psi_right[1] = psi_right[1] * matching_factor;

  double difference_in_left_and_right_derivative =
    (psi[index_classical_limit - 1] + psi_right[1] \
     - (14. - 12. * f[index_classical_limit]) * psi_matching_point_left) / delta_xi;

  printf("difference_in_left_and_right_derivative = %f\n", difference_in_left_and_right_derivative);
  
  if (fabs(difference_in_left_and_right_derivative) <= threshold_difference_in_derivatives){
    printf("difference is small enough! %f <= %f\n",
	   fabs(difference_in_left_and_right_derivative),
	   threshold_difference_in_derivatives);
    printf("returning 1\n");
    match_psi_right(matching_factor, N_intervals_right);
    return 1;
  }
  if (difference_in_left_and_right_derivative < 0) {
    epsilon_min = epsilon;
  } else {
    epsilon_max = epsilon;
  }
  return 0;
}

void match_psi_right(double matching_factor, int N_intervals_right) {
  for (int i = 2; i <= N_intervals_right; i++) {
    psi_right[i] = psi_right[i] * matching_factor;
  }  
}

void normalize_psi() {
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

void initialize_epsilon_range(void) {
  epsilon_max = vpot[N_intervals];
  epsilon_min = epsilon_max;
  
  for (int i = 0; i <= N_intervals; ++i) {
    if ( vpot[i] < epsilon_min )
      epsilon_min = vpot[i];
    if ( vpot[i] > epsilon_max )
      epsilon_max = vpot[i];
  }
}
