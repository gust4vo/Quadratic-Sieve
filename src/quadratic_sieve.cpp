#include "quadratic_sieve.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "mmc.hpp"
#include <cmath>
#include <iostream>


void quadratic_sieve(std::vector<mpz_class> &primes, mpz_class n) {

    mpz_class x = (std::sqrt(n));

    if (x*x < n) x++;

    mpz_class possible_smooth_size = primes.size() * 6;

    std::vector<mpz_class> possible_smooth(possible_smooth_size.get_ui());
    
    for (ssize_t i = 0; i < possible_smooth_size; i++)
        possible_smooth[i] = (x + i)*(x + i) % n; 
    
    std::vector<mpz_class> bases;
    bases.push_back(2);

    std::cout << "\t\t\t\t";

    for (ssize_t i = 0; i < primes.size(); i++)
    {
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
        {
            bases.push_back(primes[i]);
            std::cout << primes[i] << ' ';
        }
    }
        std::cout << '\n';

    std::vector<std::vector<unsigned long long>> pre_matriz(
    possible_smooth_size.get_ui(), 
    std::vector<unsigned long long>(bases.size(), 0));

    mpz_class tonelli_ans, ans;
    
    tonelli_ans = cipolla(n, 2);

    ans = tonelli_ans - (x % 2);

    // for (ssize_t k = 0; k < possible_smooth_size; k++)
    //     std::cout << possible_smooth[k] << ' '; 
    
    // std::cout << '\n';
    
    for (ssize_t k = ans.get_ui(); k < possible_smooth_size; k += 2)
        while(possible_smooth[k] % 2 == 0) 
        {
            possible_smooth[k] /= 2;
            pre_matriz[k][0]++;
        }

    // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[0] << '\n';

    // for (ssize_t k = 0; k < possible_smooth_size; k++)
    //     std::cout << possible_smooth[k] << ' ';

    // std::cout << '\n';
    
    for (ssize_t i = 1; i < bases.size(); i++)
    {
        tonelli_ans = cipolla(n, bases[i]);

        ans = (tonelli_ans - (x % bases[i])) % bases[i]; 

        if (ans < 0) ans += bases[i];
        
        for (ssize_t k = ans.get_ui(); k < possible_smooth_size; k += bases[i].get_ui())
        {
            while(possible_smooth[k] % bases[i] == 0)
            {
                possible_smooth[k] /= bases[i];
                pre_matriz[k][i]++;
            } 
        }

        // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[i] << '\n';
        
        // for (ssize_t k = 0; k < possible_smooth_size; k++)
        //     std::cout << possible_smooth[k] << ' ';

        // std::cout << '\n';
    
        ans = ((bases[i] - tonelli_ans) - (x % bases[i])) % bases[i];

        if (ans < 0) ans += bases[i];
        
        for (ssize_t k = ans.get_ui(); k < possible_smooth_size; k += bases[i].get_ui())
        {
            while(possible_smooth[k] % bases[i] == 0) 
            {
                possible_smooth[k] /= bases[i];
                pre_matriz[k][i]++;
            }
        }  

        // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[i] << '\n';

        // for (ssize_t k = 0; k < possible_smooth_size; k++)
        //     std::cout << possible_smooth[k] << ' ';

        // std::cout << '\n';
        
        
    }

    for (ssize_t i = 0; i < possible_smooth_size; i++)
    {
        std::cout << '\t' << (x + i)*(x + i) << "\t"  << (x + i)*(x + i) % n << "\t" << possible_smooth[i] << "\t\t";

        for (ssize_t j = 0; j < bases.size(); j++)
        {
            std::cout << pre_matriz[i][j] << ' ';
        }

        std::cout << '\n';
    }
    
        
    
    
    


}