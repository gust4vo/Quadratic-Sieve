#ifndef QUADRATIC_SIEVE
#define QUADRATIC_SIEVE

#include <gmpxx.h>
#include <vector>

void quadratic_sieve(const std::vector<mpz_class> &primes, const mpz_class &n);


#endif