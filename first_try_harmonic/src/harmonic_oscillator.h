#ifndef HARMONIC_H
#define HARMONIC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double omega;
double epsilon;
double delta_xi;
double delta_xi_squared;

double *psi, *xi, *f, *vpot;

double V(double xi);

double K(void);

double psi_n_plus(double xi_n_minus, double psi_n_minus,
		  double xi_n, double psi_n,
		  double xi_n_plus);

double get_g_n(double xi_n);

double get_f_n(double g_n);

double get_psi_1_for_nodes_even(double psi_0);

int compute_psi(double epsilon_min, double epsilon_max,
		 int N_intervals,
		 int n_max_iterations,
		 int expected_number_of_nodes);

void compute_f(int N_intervals);

#endif /* HARMONIC_H */
