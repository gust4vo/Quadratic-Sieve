#ifndef ERATHOSTENES_HPP
#define ERATHOSTENES_HPP

#include <vector>
#include <gmpxx.h>

void GetPrimes(std::vector<mpz_class>& primes, unsigned long int upperBound);

#endif