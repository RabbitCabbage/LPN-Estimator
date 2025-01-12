#include <stdio.h>

#include <gmp.h>
#include <mpfr.h>

#define PRECISION 200

// Combination
void combination(mpfr_t &C, unsigned int n, unsigned int k)
{

  unsigned int i;
  mpfr_init2(C, PRECISION);
  if (k > n) {
    mpfr_set_ui(C, 0, MPFR_RNDN);
    return;
  }
  mpfr_set_ui(C, 1, MPFR_RNDN);

  if (k > n / 2)
    k = n - k;

  for (i = 1; i <= k; i++)
  {
    mpfr_mul_ui(C, C, n - k + i, MPFR_RNDN);
    mpfr_div_ui(C, C, i, MPFR_RNDN);
  }
}

// Pooled Gaussian, cost in log2
void Gauss(mpfr_t &T, unsigned int N, unsigned int k, unsigned int t)
{
  // prepare variables, strassen cst = 2.8
  mpfr_t total_, success_, trials_, strassen_, cst_;
  mpfr_init2(total_, PRECISION);
  mpfr_init2(success_, PRECISION);
  mpfr_init2(trials_, PRECISION);
  mpfr_init2(strassen_, PRECISION);
  mpfr_init2(cst_, PRECISION);
  mpfr_set_d(cst_, 2.8, MPFR_RNDN);

  // expected trials in log2
  combination(total_, N, k);
  combination(success_, N - k, t);
  mpfr_div(trials_, total_, success_, MPFR_RNDN);
  mpfr_log2(trials_, trials_, MPFR_RNDN);

  // cost of Strassen's algorithm
  mpfr_pow(strassen_, ((N - k > k) ? k : N - k), cst_, MPFR_RNDN);
  mpfr_log2(strassen_, strassen_, MPFR_RNDN);

  // total cost
  mpfr_add(T, trials_, strassen_, MPFR_RNDN);

  // clear variables
  mpfr_clear(total_);
  mpfr_clear(trials_);
  mpfr_clear(strassen_);
  mpfr_clear(cst_);
  mpfr_free_cache();
}

int main()
{
  unsigned int i;
  mpfr_t s, t, u;

  mpfr_init2(t, 200);
  mpfr_set_d(t, 1.0, MPFR_RNDN);
  mpfr_init2(s, 200);
  mpfr_set_d(s, 1.0, MPFR_RNDN);
  mpfr_init2(u, 200);
  for (i = 1; i <= 100; i++)
  {
    mpfr_mul_ui(t, t, i, MPFR_RNDN);
    mpfr_set_d(u, 1.0, MPFR_RNDN);
    mpfr_div(u, u, t, MPFR_RNDN);
    mpfr_add(s, s, u, MPFR_RNDN);
  }
  printf("Sum is ");
  mpfr_out_str(stdout, 10, 0, s, MPFR_RNDN);
  putchar('\n');
  mpfr_clear(s);
  mpfr_clear(t);
  mpfr_clear(u);
  mpfr_free_cache();
  return 0;
}
