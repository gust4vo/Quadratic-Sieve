#include "cipolla.hpp"
#include "ak_mod_n.hpp"

Complex complexMul(const Complex& x, const Complex& y, const mpz_class& n, const mpz_class& w) {
    return Complex((x.a * y.a + x.b * y.b * w) % n, (x.a * y.b + x.b * y.a) % n);
}

Complex complexPow(Complex base, mpz_class exp, const mpz_class& n, const mpz_class& w) {
    Complex result(1, 0);
    while (exp > 0) {
        if (exp % 2 == 1) {
            result = complexMul(result, base, n, w);
        }
        base = complexMul(base, base, n, w);
        exp /= 2;
    }
    return result;
}

mpz_class cipolla(mpz_class n, mpz_class p) {

    if (p == 2) return n % 2;
    
    if (n == 0) return 0;
    if (mpz_legendre(n.get_mpz_t(), p.get_mpz_t()) != 1) {
        throw std::invalid_argument("No solution exists");
    }

    mpz_class a = 0, w;
    while (true) {
        a++;
        w = (a * a - n) % p;
        if (mpz_legendre(w.get_mpz_t(), p.get_mpz_t()) == -1) break;
    }

    Complex x(a, 1);
    mpz_class exp = (p + 1) / 2;
    Complex res = complexPow(x, exp, p, w);

    if (res.a < 0) return (p - res.a) % p;
    return res.a;
}