#ifndef TONELLI_SHANKS_HPP
#define TONELLI_SHANKS_HPP

#include <gmpxx.h>

mpz_class choose_non_residue(const mpz_class p);

mpz_class tonelli_shanks(const mpz_class& a, const mpz_class& p);

mpz_class tonelli_shanks_recursive(const mpz_class& a, const mpz_class& b, const mpz_class& m, const mpz_class& p);

mpz_class tonelli_shanks_iterative(const mpz_class& a, const mpz_class& p);

#endif