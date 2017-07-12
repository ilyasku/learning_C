#ifndef HARMONIC_H
#define HARMONIC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

double sqrt();

double omega;
double epsilon;
double delta_xi;
double delta_xi_squared;
double xi_classical_limit;

double epsilon_min, epsilon_max;

double *psi, *xi, *f, *vpot;
double *psi_right;

double psi_matching_point_left, psi_matching_point_right;

double V(double xi);

double K(void);

double psi_n_plus(double xi_n_minus, double psi_n_minus,
		  double xi_n, double psi_n,
		  double xi_n_plus);

double get_g_n(double xi_n);

double get_f_n(double g_n);

double get_psi_1_for_nodes_even(double psi_0);

int compute_psi(int N_intervals,
		int n_max_iterations,
		int expected_number_of_nodes);

void compute_f(int N_intervals);

int find_classical_limit(int lower_limit, int upper_limit);

int compute_psi_from_right_to_left(int N_intervals, int index_classical_limit,
				   int n_max_iterations,
				   double threshold_difference_in_derivatives);

void normalize_psi_right_and_write_to_psi(double matching_factor, int N_intervals_right,
					  int index_classical_limit);

void normalize_psi(int N_intervals);

#endif /* HARMONIC_H */
