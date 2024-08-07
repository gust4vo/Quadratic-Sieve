#include "cipolla.hpp"
#include "ak_mod_n.hpp"

mpz_class cipolla_iterative(const mpz_class& n, const mpz_class& p) {
    if (ak_mod_n(n, (p - 1) / 2, p) != 1) {
        return -1; // n não é um resíduo quadrático
    }

    // Encontrar a e ω tal que (a^2 - n) é um não-resíduo quadrático
    mpz_class a, omega, w;
    do {
        a = rand() % p; // Escolher um valor aleatório para a
        w = (a * a - n) % p;
        omega = (p - w) % p;
    } while (ak_mod_n(omega, (p - 1) / 2, p) == 1);

    // Inicializar as variáveis para a parte iterativa
    mpz_class x = 1, y = 0, v = a, u = 1;
    mpz_class m = (p + 1) / 2;

    while (m > 0) {
        if (m % 2 == 1) {
            mpz_class temp_x = (x * v + y * u * omega) % p;
            y = (x * u + y * v) % p;
            x = temp_x;
        }
        mpz_class temp_v = (v * v + u * u * omega) % p;
        u = (2 * v * u) % p;
        v = temp_v;
        m /= 2;
    }

    return x;
}

Cipolla cipolla_recursive(const Cipolla& x, mpz_class n) {
    if (n == 0) return Cipolla(1, 0, x.omega, x.p);
    Cipolla half = cipolla_recursive(x, n / 2);
    half = half * half;
    if (n % 2 == 1) half = half * x;
    return half;
}

mpz_class cipolla(const mpz_class& n, const mpz_class& p_input) {
    mpz_class p = p_input;

    if (ak_mod_n(n, (p - 1) / 2, p) != 1) {
        return -1; // n não é um resíduo quadrático
    }

    mpz_class a, omega, w;
    do {
        a = rand() % p;
        w = (a * a - n) % p;
        omega = (p - w) % p;
    } while (ak_mod_n(omega, (p - 1) / 2, p) == 1);

    Cipolla x(a, 1, omega, p);
    Cipolla result = cipolla_recursive(x, (p + 1) / 2);

    return result.a;
}