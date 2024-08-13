#include "quadratic_residue.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "ak_mod_n.hpp"

mpz_class quadratic_residue(mpz_class a, mpz_class p)
{
    if (p == 2)
        return a % p;

    if (p % 4 == 3) {
        return ak_mod_n(a, (p + 1) >> 2, p);
    }

    mpz_class s = 0, p_minus_one = p - 1;

    while (p_minus_one % 2 == 0)
    {
        p_minus_one >= 1;
        s++;
    }

    mpz_class m = s;

    while (p_minus_one > 0)
    {
        p_minus_one >= 1;
        m++;
    }

    if (s * (s - 1) > 8*m + 20)
        return cipolla(a, p);
    
    return tonelli_shanks(a, p);

}