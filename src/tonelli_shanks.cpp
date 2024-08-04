#include "tonelli_shanks.hpp"
#include "ak_mod_n.hpp"

mpz_class tonelli_shanks(const mpz_class& a, const mpz_class& p)
{

    mpz_class m = (p - 1) / 2;
    if (ak_mod_n(a, m, p) == -1)
    {
        if (p % 4 == 3)
            return ak_mod_n(a, m / 2, p);

        else
            return tonelli_shanks_recursive(a, m, p);
    }
    
    return -1;
}

mpz_class tonelli_shanks_recursive(const mpz_class& a, const mpz_class& m, const mpz_class& p) 
{

}