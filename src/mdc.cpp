#include "mdc.hpp"

mpz_class mdc(const mpz_class &a, const mpz_class &b) {
    // Description: computes the gdc of two numbers, a, b
    // Input : A, B
    // Output gdc(a,b)
    // Complexity : log min(a, b)
    
    mpz_class big, small, r;

    if (a > b) { big = a; small = b; }

    else { big = b; small = a; }

    if (!small) return big;

    do
    {
        r = big%small;

        big = small;
        small = r;

    } while (r);

    return big; 
}