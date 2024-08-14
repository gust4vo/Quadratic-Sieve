#include <iostream>
#include <ctime>
#include <iomanip>
#include <cmath>
#include <vector>
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "quadratic_sieve.hpp"
#include "erathostenes.hpp"
#include "gauss_jordan.hpp"

void calculateTimeDifference(struct timespec t1, struct timespec t2, struct timespec * res)
{
    if (t2.tv_nsec < t1.tv_nsec) {
        res->tv_nsec = 1000000000 + t2.tv_nsec - t1.tv_nsec;
        res->tv_sec = t2.tv_sec - t1.tv_sec - 1;
    } else {
        res->tv_nsec = t2.tv_nsec - t1.tv_nsec;
        res->tv_sec = t2.tv_sec - t1.tv_sec;
    }
}


int main()
{
    std::cout << "===============================================================================================================\n";
    std::cout << "                                                 CRIVO QUADRÁTICO                                              \n";
    std::cout << "===============================================================================================================\n\n";

    mpz_class n;
    std::vector<mpz_class> primes;
    std::cout << "N = "; std::cin >> n;

    unsigned long int upperBound = exp(0.58 * sqrt(log(n.get_d()) * log(log(n.get_d()))));

    std::cout << "\nLimitante superior (B) = " << upperBound << std::endl;

    GetPrimes(primes, upperBound);

    std::cout << "\nQuantidade de primos encontrados : " << primes.size() << std::endl;

    struct timespec initTime, endTime, diffTime;
    clock_gettime(CLOCK_MONOTONIC, &initTime);
    quadratic_sieve(primes, n);
    clock_gettime(CLOCK_MONOTONIC, &endTime);

    calculateTimeDifference(initTime, endTime, &diffTime);

    std::cout << "\nTempo levado: " << diffTime.tv_sec << "." << std::setw(9) << std::setfill('0') << diffTime.tv_nsec << "s\n\n";

    return 0;
}