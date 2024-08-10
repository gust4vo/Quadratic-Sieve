#include "erathostenes.hpp"

void Erathostenes::GetPrimes(std::vector<unsigned long int>& primes, unsigned long int upperBound) {
    primes_bitset.resize(upperBound + 1, true); 
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
