#ifndef CIPOLLA_HPP
#define CIPOLLA_HPP

#include <gmpxx.h>

struct Complex {
    mpz_class a, b;
    Complex(mpz_class a = 0, mpz_class b = 0) : a(a), b(b) {}
};

Complex complexMul(const Complex& x, const Complex& y, const mpz_class& n, const mpz_class& w);

Complex complexPow(Complex base, mpz_class exp, const mpz_class& n, const mpz_class& w);

mpz_class cipolla(const mpz_class &n, const mpz_class &p);

#endif
