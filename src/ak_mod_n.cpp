#include "ak_mod_n.hpp"

mpz_class ak_mod_n(const mpz_class& a, const mpz_class& k, const mpz_class& n) {

    // Description : computes a powered k mod n
    // Input : a, k, n
    // Output : module a powered k of n
    // Complexity : O(log n)

    mpz_class number, pow, ans;
    number =  a % n;
    ans = 1; pow = k;

    if ((pow & 1) == 1) ans = number;
    pow >>= 1;

    while (pow != 0)
    {
        number = (number*number) % n;

        if ((pow & 1) == 1) ans = (number*ans) % n;
    
        pow >>= 1;
    }

    return ans;

}