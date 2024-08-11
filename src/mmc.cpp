#include "mmc.hpp"
#include "mdc.hpp"

mpz_class mmc(mpz_class a, mpz_class b) {
    // Description : computes the lsn of two numbers, a, b
    // Input : A, B
    // Output : lsn(A.B)
    // Complexity : log min(a, b)
    return (a*b) / mdc(a, b);
}