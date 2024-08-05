#include "tonelli_shanks.hpp"
#include "ak_mod_n.hpp"
#include "inverse_a_modulo_n.hpp"

mpz_class choose_non_residue(const mpz_class p) {

    mpz_class m = (p - 1) / 2;
    for (mpz_class i = 2; i < p; i++)
    {
        if (ak_mod_n(i, m, p) == p - 1)
            return i;
    }
    
}

mpz_class tonelli_shanks(const mpz_class& a, const mpz_class& p)
{

    mpz_class m = (p - 1) >> 1;

    if (ak_mod_n(a, m, p) == 1)
    {
        if (p % 4 == 3)
            return ak_mod_n(a, (p + 1) >> 2, p);

        else
            return tonelli_shanks_recursive(a, choose_non_residue(p) , m, p) % p;
    }
    
    return -1;
}

mpz_class tonelli_shanks_recursive(const mpz_class& a, const mpz_class& b, const mpz_class& m, const mpz_class& p) 
{
    if (ak_mod_n(a, m, p) == p - 1)
        return  (inverse_a_modulo_n(b % p, p) * tonelli_shanks_recursive(a*b*b % p, b*b % p, m, p)) % p;

    else {

        if (m % 2 == 1)
            return ak_mod_n(a, (m + 1) / 2, p);
        
        else
            return tonelli_shanks_recursive(a, b, m / 2, p);

    }
}   


mpz_class tonelli_shanks_iterative(const mpz_class& a, const mpz_class& p)
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