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

  N_intervals = atoi(argv[1]);
  delta_xi = atof(argv[2]); 
  double psi_one_for_nodes_odd = atof(argv[3]);
  double psi_zero_for_nodes_even = atof(argv[4]);  
  char *output_folder_name = argv[5];

  int n_max_iterations = 1000;
  double threshold_in_derivatives = 1e-5;

  psi = (double *) malloc((N_intervals + 1) * sizeof(double));
  xi = (double *) malloc((N_intervals + 1) * sizeof(double));
  f = (double *) malloc((N_intervals + 1) * sizeof(double));
  vpot = (double *) malloc((N_intervals + 1) * sizeof(double));

  PSI_RIGHT_ALLOCATED = 0;
  
  /* set up the potential */
  for (int i = 0; i <= N_intervals; i++) {
    xi[i] = (double) i * delta_xi;
    vpot[i] = 0.5 * xi[i] * xi[i];
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

    printf("======================================\n");
    printf("number_of_nodes = %i\n", number_of_nodes);
    printf("psi_0 = %f\n", psi_0);
    printf("psi_1 = %f\n", psi_1);

    
    initialize_epsilon_range();
    
    psi[0] = psi_0;
    psi[1] = psi_1;
            
    int remaining_iterations = compute_psi(n_max_iterations, number_of_nodes, threshold_in_derivatives);

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

    fprintf(stderr, "------------------\n");
    fprintf(stderr, "before writing the file:\n");
    fprintf(stderr, "N_intervals = %i\n", N_intervals);	
    
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
