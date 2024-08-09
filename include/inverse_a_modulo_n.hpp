#ifndef INVERSE_A_MODULO_N
#define INVERSE_A_MODULO_N

#include <vector>
#include <gmpxx.h>

mpz_class inverse_a_modulo_n(const mpz_class& a, const mpz_class& n);

#endif