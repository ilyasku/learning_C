#ifndef HARMONIC_H
#define HARMONIC_H

double K;
double epsilon;
double delta_xi_squared;

double V(double x);

double psi_n_plus(double xi_n_minus, double psi_n_minus,
		  double xi_n, double psi_n,
		  double xi_n_plus);

#endif /* HARMONIC_H */
