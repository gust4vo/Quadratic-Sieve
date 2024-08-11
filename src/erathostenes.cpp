#include "erathostenes.hpp"

void GetPrimes(std::vector<mpz_class>& primes, unsigned long int upperBound) {
    std::vector<bool> primes_bitset(upperBound + 1, true);
    primes.resize(upperBound);

    primes_bitset[0] = primes_bitset[1] = false; 

    for (unsigned long int p = 2; p * p <= upperBound; ++p) {
        if (primes_bitset[p]) {
            for (unsigned long int i = p * p; i <= upperBound; i += p)
                primes_bitset[i] = false;
        }
    }

    for (unsigned long int p = 2; p <= upperBound; ++p) {
        if (primes_bitset[p]) 
            primes.push_back(p);

    }
}
