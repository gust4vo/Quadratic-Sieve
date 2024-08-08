#include "tonelli_shanks.hpp"
#include "ak_mod_n.hpp"

mpz_class tonelli_shanks(const mpz_class& a, const mpz_class& p)
{

    mpz_class k = (p - 1) >> 1;

    if(ak_mod_n(a, k, p) != 1) {
        return -1;
    }

    if (p % 4 == 3) {
        return ak_mod_n(a, k >> 1, p);
    }

    mpz_class s = 0;
    mpz_class q = p-1;

    // encontrar expoente
    while (q % 2 == 0) {
        q >>= 1;
        s++;
    }

    // encontrar z
    mpz_class z = 2;
    while (ak_mod_n(z, k, p) == 1) {
        z++;
    }

    mpz_class m = s;
    mpz_class c = ak_mod_n(z, q, p);
    mpz_class t = ak_mod_n(a, q, p);
    mpz_class r = ak_mod_n(a, (q + 1) >> 1, p);

    while (t != 1) {
        mpz_class t2 = t;
        mpz_class i = 0;
        while (t2 != 1) {
            t2 = ak_mod_n(t2, 2, p);
            i++;
        }

        mpz_class b = ak_mod_n(c, ak_mod_n(2, m - i - 1, p), p);
        m = i;
        c = ak_mod_n(b, 2, p);
        t = (t * c) % p;
        r = (r * b) % p;
    }

    return r;

}