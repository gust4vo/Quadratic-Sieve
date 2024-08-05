#include "inverse_a_modulo_n.hpp"
#include "euclidean_extension.hpp"

mpz_class inverse_a_modulo_n(const mpz_class& a, const mpz_class& n)
{
    // Description : computes the inverse of a modulo n
    // Input : a, n
    // Output : inverse of a modulo n
    // Complexity: O(log(min(a, n))

    if(a == 1)
        return 1;

    mpz_class result;
    if(a > n) result = euclidean_extension(a, n).first;
    else result = euclidean_extension(n, a).second;

    result = (result % n + n) % n;

    return result;
}