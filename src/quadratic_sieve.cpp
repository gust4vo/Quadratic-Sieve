#include "quadratic_sieve.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "mdc.hpp"
#include "ak_mod_n.hpp"
#include "gauss_jordan.hpp"
#include <cmath>
#include <iostream>

void sieving(std::vector<mpz_class> possible_smooth) {

}

void quadratic_sieve(std::vector<mpz_class> &primes, mpz_class n) {
    mpz_class x = (std::sqrt(n));

    if (x*x < n) x++;   

    size_t possible_smooth_size = primes.size()*6;

    std::vector<mpz_class> possible_smooth(possible_smooth_size);

    for (size_t i = 0; i < possible_smooth_size; i++) {
        possible_smooth[i] = (x + i) * (x + i) % n;
    }
    
    std::vector<mpz_class> bases;
    bases.push_back(2);

    for (size_t i = 0; i < primes.size(); i++)
    {
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
        {
            bases.push_back(primes[i]);
        }
    }

    std::vector<std::vector<unsigned long long>> pre_matriz(
    possible_smooth_size,
    std::vector<unsigned long long>(bases.size(), 0));

    mpz_class tonelli_ans, ans;
    
    tonelli_ans = cipolla(n, 2);

    ans = tonelli_ans - (x % 2);
    for (size_t k = ans.get_ui(); k < possible_smooth_size; k += 2) {
        while(possible_smooth[k] % 2 == 0) 
        {
            possible_smooth[k] /= 2;
            pre_matriz[k][0]++;
        }
    }

    for (size_t i = 1; i < bases.size(); i++)
    {
        tonelli_ans = cipolla(n, bases[i]);

        ans = (tonelli_ans - (x % bases[i])) % bases[i]; 

        if (ans < 0) ans += bases[i];
        
        for (size_t k = ans.get_ui(); k < possible_smooth_size; k += bases[i].get_ui())
        {
            while(possible_smooth[k] % bases[i] == 0)
            {
                possible_smooth[k] /= bases[i];
                pre_matriz[k][i]++;
            } 
        }
    
        ans = ((bases[i] - tonelli_ans) - (x % bases[i])) % bases[i];

        if (ans < 0) ans += bases[i];
        
        for (size_t k = ans.get_ui(); k < possible_smooth_size; k += bases[i].get_ui())
        {
            while(possible_smooth[k] % bases[i] == 0) 
            {
                possible_smooth[k] /= bases[i];
                pre_matriz[k][i]++;
            }
        }  
    }

    std::vector<size_t> smooth_index;
    std::vector<std::vector<unsigned long long>> linear_system;
    
    for (size_t i = 0, aux = 0; i < possible_smooth_size; i++) {
        if (possible_smooth[i] == 1 && aux < bases.size() + 1) {
            smooth_index.push_back(i);
            linear_system.push_back(pre_matriz[i]);
            aux++;
        }

        else
        {
            pre_matriz[i].clear();
            pre_matriz[i].shrink_to_fit();
        }
    }

    possible_smooth.clear();
    possible_smooth.shrink_to_fit();

    if (smooth_index.size() < bases.size() + 1) {
        std::cout << "There's no sufficient smooth numbers" << std::endl;
        return;
    }
    
    // resolve the linear system                           
    std::vector<std::vector<int>> solutions = gauss_jordan(linear_system);

    std::cout << "Oi\n";
    mpz_class a = 1, b = 1;
    for(size_t combination=0; combination < solutions.size(); combination++) {
        for (size_t i = 0; i < smooth_index.size(); i++) {
            if (solutions[combination][i]) {
                for(size_t j = 0; j < bases.size(); j++) {
                    mpz_class p;
                    mpz_pow_ui(p.get_mpz_t(), bases[j].get_mpz_t(), pre_matriz[smooth_index[i]][j]);
                    a *= p;
                }
                b *= x + smooth_index[i];
            }
        }

        mpz_sqrt(a.get_mpz_t(), a.get_mpz_t());

        if(a%n != b%n) break;

        a = 1, b = 1;
    }

    mpz_class f1 = mdc(b - a, n);
    mpz_class f2 = mdc(a + b, n);
    std::cout << "n = " << f1 << " * " << f2 << "\n";
}