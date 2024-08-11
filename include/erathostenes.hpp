#ifndef ERATHOSTENES_HPP
#define ERATHOSTENES_HPP

#include <vector>

class Erathostenes {
private:
    std::vector<bool> primes_bitset;

public:
    void GetPrimes(std::vector<unsigned long int>& primes, unsigned long int upperBound);
};


#endif