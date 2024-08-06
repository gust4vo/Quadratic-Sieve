#ifndef EUCLIDEAN_EXTENSION
#define EUCLIDEAN_EXTENSION

#include <gmpxx.h>

std::pair<mpz_class, mpz_class> euclidean_extension(const mpz_class& a, const mpz_class& b);

#endif