#include "quadratic_sieve.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include <cmath>
#include <iostream>


void quadratic_sieve(std::vector<mpz_class> &primes, mpz_class n) {

    mpz_class x = std::sqrt(n);

    std::vector<mpz_class> possible_smooth(primes.size()*6);
    
    for (ssize_t i = 0; i < primes.size()*6; i++)
        possible_smooth[i] = (x + i)*(x + i) % n; 
    

    std::vector<mpz_class> bases;
    bases.push_back(2);

    for (ssize_t i = 0; i < primes.size(); i++)
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
            bases.push_back(primes[i]);
    
    for (ssize_t i = 1; i < bases.size(); i++)
    {
        for (ssize_t j = 0; j < bases[i]; j++)
        {
            mpz_class tonelli_ans = cipolla(n, bases[i]);

            mpz_class ans = tonelli_ans - (x % bases[i]); 
            
            for (ssize_t k = 0; k < primes.size()*6; k += ans.get_ui())
            {
                if (possible_smooth[k] != 1 && bases[i] <= possible_smooth[k]) possible_smooth[k] /= bases[i];
            }
            ans = (bases[i] - tonelli_ans) - (x % bases[i]);
            
            for (ssize_t k = 0; k < primes.size()*6; k += ans.get_ui())
                if (possible_smooth[k] != 1 && bases[i] <= possible_smooth[k]) possible_smooth[k] /= bases[i];
        }
        
    }




    for (ssize_t i = 0; i < primes.size()*3; i++)
        std::cout << (x + i)*(x + i) % n << " " << possible_smooth[i*2] << '\n';
    
    
    


}