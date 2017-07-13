#ifndef HARMONIC_H
#define HARMONIC_H

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

double sqrt();

int N_intervals;
int index_classical_limit;

int PSI_RIGHT_ALLOCATED;

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

int compute_psi(int n_max_iterations,
		int expected_number_of_nodes,
		double threshold_difference_in_derivatives);

void compute_f();

void concatenate_psi_left_and_psi_right(void);

void find_classical_limit();

int compute_psi_right(double threshold_difference_in_derivatives);

void match_psi_right(double matching_factor, int N_intervals_right);
			

void normalize_psi();

void initialize_epsilon_range(void);

#endif /* HARMONIC_H */
