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