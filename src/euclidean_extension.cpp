#include "euclidean_extension.hpp"

std::pair<mpz_class, mpz_class> euclidean_extension(const mpz_class& a, const mpz_class& b) 
{
    // Description : Computes the minimal solution of the diophantine equation a*x + b*y = 1
    // Input : a, b
    // Output : x, y
    // Complexity : O(log (min(a, b))

    mpz_class big, small, x, x0, x1, y, y0, y1, q, r;

    big = a; small = b;
    x = 1; x0 = 1; x1 = 0;
    y = 0; y0 = 0; y1 = 1;

    while (true)
    {   
        q = big / small;
        r = big - (q*small);
        
        if (r == 0) break;

        x = x0 - x1 * q;
        y = y0 - y1 * q;

        x0 = x1; x1 = x;
        y0 = y1; y1 = y;

        big = small;
        small = r;
    }
    
    return std::pair<mpz_class, mpz_class> {x, y};
}