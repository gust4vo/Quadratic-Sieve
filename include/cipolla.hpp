#ifndef CIPOLLA_HPP
#define CIPOLLA_HPP

#include <iostream>
#include <gmpxx.h>

struct Cipolla {
    mpz_class a, b, omega, p;
    Cipolla(const mpz_class& _a, const mpz_class& _b, const mpz_class& _omega, const mpz_class& _p)
        : a(_a), b(_b), omega(_omega), p(_p) {}

    Cipolla operator*(const Cipolla& other) const {
        mpz_class new_a = (a * other.a + b * other.b * omega) % p;
        mpz_class new_b = (a * other.b + b * other.a) % p;
        return Cipolla(new_a, new_b, omega, p);
    }
};

mpz_class cipolla_iterative(const mpz_class& n, const mpz_class& p);
Cipolla cipolla_recursive(const Cipolla& x, mpz_class n);
mpz_class cipolla(const mpz_class& n, const mpz_class& p_input);

#endif
