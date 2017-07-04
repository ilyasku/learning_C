#include <stdlib.h>
#include <check.h>
#include "src/harmonic_oscillator.h"

START_TEST(test_psi_n_plus)
{
  epsilon = 1.0;
  delta_xi = 0.1;
  delta_xi_squared = delta_xi * delta_xi;

  double psi_0 = 1.0;
  double psi_1 = 2.0;

  double xi_0 = 0.0;
  double xi_1 = delta_xi;
  double xi_2 = delta_xi * 2;

  double psi_2 = psi_n_plus(xi_0, psi_0, xi_1, psi_1, xi_2);
  ck_assert_double_eq_tol(3.039798, psi_2, 0.00001);
}
END_TEST

START_TEST(test_get_psi_1_for_nodes_odd)
{
  epsilon = 1.0;
  delta_xi = 0.1;
  delta_xi_squared = delta_xi * delta_xi;

  double psi_0 = 1.0;

  double psi_1 = get_psi_1_for_nodes_odd(psi_0);
  ck_assert_double_eq_tol(1.010008, psi_1, 0.00001);
}
END_TEST

Suite * harmonic_suite(void)
{
  Suite *s;
  TCase *tc_core;
  s = suite_create("Harmonic");

  tc_core = tcase_create("Core");
  
  tcase_add_test(tc_core, test_psi_n_plus);
  tcase_add_test(tc_core, test_get_psi_1_for_nodes_odd);
  suite_add_tcase(s, tc_core);

  return s;
}

int main(void) {
  int number_failed;
  Suite *s;
  SRunner *sr;

  s = harmonic_suite();
  sr = srunner_create(s);

  srunner_run_all(sr, CK_NORMAL);
  number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);
  return (number_failed == 0) ? EXIT_SUCCESS: EXIT_FAILURE;
}

