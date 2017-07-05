#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "harmonic_oscillator.h"


int main(int argc, char **argv) {  
  if (argc != 6)
  {
    fprintf(stderr, "need 5 arguments!\n");
    abort();
  } 
  
  FILE *out_file;

  int N_intervals = atoi(argv[1]);
  delta_xi = atof(argv[2]); 
  double psi_one_for_nodes_odd = atof(argv[3]);
  double psi_zero_for_nodes_even = atof(argv[4]);  
  char *output_folder_name = argv[5];

  int n_max_iterations = 1000;

  psi = (double *) malloc((N_intervals + 1) * sizeof(double));
  xi = (double *) malloc((N_intervals + 1) * sizeof(double));
  f = (double *) malloc((N_intervals + 1) * sizeof(double));
  vpot = (double *) malloc((N_intervals + 1) * sizeof(double));

  /* set up the potential */
  for (int i = 0; i <= N_intervals; i++) {
    xi[i] = (double) i * delta_xi;
    vpot[i] = 0.5 * xi[i] * xi[i];
  }  
  
  double epsilon_max = vpot[N_intervals];
  double epsilon_min = epsilon_max;

  for (int i = 0; i <= N_intervals; ++i) {
    if ( vpot[i] < epsilon_min )
      epsilon_min = vpot[i];
    if ( vpot[i] > epsilon_max )
      epsilon_max = vpot[i];
  }
  
  double psi_0;
  double psi_1;

  double sign_for_negative_xi;
  
  delta_xi_squared = delta_xi * delta_xi;
  
  for (int number_of_nodes = 0; number_of_nodes < 4; number_of_nodes ++) {    
    if ((number_of_nodes % 2) == 0) {
      psi_0 = psi_zero_for_nodes_even;
      psi_1 = get_psi_1_for_nodes_even(psi_0);
      sign_for_negative_xi = 1.;
    } else {
      psi_0 = 0.0;
      psi_1 = psi_one_for_nodes_odd;
      sign_for_negative_xi = -1.;
    }
    printf("number_of_nodes = %i\n", number_of_nodes);
    printf("psi_0 = %f\n", psi_0);
    printf("psi_1 = %f\n", psi_1);

    psi[0] = psi_0;
    psi[1] = psi_1;

    int remaining_iterations = compute_psi(epsilon_min, epsilon_max, N_intervals,
					   n_max_iterations, number_of_nodes);

    printf("number of actual iterations: %d\n", n_max_iterations - remaining_iterations);
    
    char *file_name;
    file_name = malloc(strlen(output_folder_name) + 1 + 6);
    strcpy(file_name, output_folder_name);
    strcat(file_name, "/out_");
    char index[2];
    sprintf(index, "%d", number_of_nodes);
    strcat(file_name, index);

    out_file = fopen(file_name, "w");

    fprintf (out_file,"#   xi       Psi(xi)      V\n");

    for (int i = N_intervals; i >=0; i--) {
      fprintf(out_file, "%7.3f%16.8e%12.6f\n",
	      -xi[i], sign_for_negative_xi* psi[i], vpot[i]);
    }
    
    for (int i = 0; i <= N_intervals; ++i) {
      fprintf(out_file, "%7.3f%16.8e%12.6f\n",
	      xi[i], psi[i], vpot[i]);
    }
    
    fclose(out_file);
    free(file_name);				  
  }  
  
  free(psi);
  free(xi);
  free(f);
  free(vpot);  
}
