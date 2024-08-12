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

#ifdef TEST_TONELLI_CIPOLLAA

int main()
{
    mpz_class n, p;
    std::cin >> n >> p;

    struct timespec initTime, endTime, diffTime;

    // Tonelli Shanks tests

    clock_gettime(CLOCK_MONOTONIC, &initTime);
    mpz_class ans2 = tonelli_shanks(n, p);
    clock_gettime(CLOCK_MONOTONIC, &endTime);

    calculateTimeDifference(initTime, endTime, &diffTime);

    std::cout << "Tonelli-Shanks: "  << ans2 << '\n';
    std::cout << "Tempo levado: " << diffTime.tv_sec << "." << std::setw(9) << std::setfill('0') << diffTime.tv_nsec << "s" << std::endl;

    // Cipolla tests
    
    clock_gettime(CLOCK_MONOTONIC, &initTime);
    ans2 = cipolla(n, p);
    clock_gettime(CLOCK_MONOTONIC, &endTime);

    calculateTimeDifference(initTime, endTime, &diffTime);
    std::cout << "Cipolla: " << ans2 << '\n';
    std::cout << "Tempo levado: " << diffTime.tv_sec << "." << std::setw(9) << std::setfill('0') << diffTime.tv_nsec << "s" << std::endl;
}

#endif

#ifdef QUADRATIC_SIEVE

int main()
{

    
    mpz_class n;
    std::vector<mpz_class> primes;
    std::cin >> n;

    unsigned long int upperBound = exp(1 * sqrt(log(n.get_d())*log(log(n.get_d()))));

    GetPrimes(primes, upperBound);

    struct timespec initTime, endTime, diffTime;

    clock_gettime(CLOCK_MONOTONIC, &initTime);
    quadratic_sieve(primes, n);
    clock_gettime(CLOCK_MONOTONIC, &endTime);

    calculateTimeDifference(initTime, endTime, &diffTime);

    
    std::cout << "Tempo levado: " << diffTime.tv_sec << "." << std::setw(9) << std::setfill('0') << diffTime.tv_nsec << "s" << std::endl;

//     std::vector<std::vector<unsigned long long>> matriz = {{0, 1, 1, 0}, {1, 0, 1, 0}, {1, 1, 0, 1}, {1, 1, 1, 0}, {0, 1, 0, 1}};
//     std::vector<int> sol = gauss_jordan(matriz);
//     // for(int i : sol) std::cout << i << " ";
//     // std::cout << "\n";
}

#endif