#include "quadratic_sieve.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "mdc.hpp"
#include "ak_mod_n.hpp"
#include "gauss_jordan.hpp"
#include <cmath>
#include <iostream>
#include <omp.h>
#include <utility>


void sieving(std::vector<mpz_class> &bases, mpz_class &n, std::vector<std::pair<mpz_class, mpz_class>> &smooth_index, std::vector<std::vector<unsigned long long>> &matrix) {

    mpz_class x = (std::sqrt(n));

    if (x*x < n) x++;

    mpz_class tonelli_ans, ans;
    
    tonelli_ans = cipolla(n, 2);
    ans = tonelli_ans - (x % 2);

    std::vector<std::pair<mpz_class, mpz_class>> possible_smooth(bases.size());
    possible_smooth[0] = {ans, ans};

    for (size_t i = 1; i < bases.size(); i++)
    {
        tonelli_ans = cipolla(n, bases[i]);

        ans = (tonelli_ans - (x % bases[i])) % bases[i]; 

        if (ans < 0) ans += bases[i];
        
        possible_smooth[i].first = ans;
    
        ans = ((bases[i] - tonelli_ans) - (x % bases[i])) % bases[i];

        if (ans < 0) ans += bases[i];
        
        possible_smooth[i].second = ans;
    }
    

    
    // omp_set_num_threads(omp_get_max_threads());  // Ajuste o número de threads, por exemplo, para 4

    bool done = false;
    
   
    

     size_t offset = 0;

    while(!done)
    {

        #pragma omp parallel
        {
            
            #pragma omp for
            for (size_t i = 0; i < bases.size(); i++)
            {
                if (!done)
                {
                    std::pair<mpz_class, mpz_class> smooth((x + i + offset), (x + i + offset) * (x + i + offset) % n);
                    mpz_class smooth_residue = smooth.second;

                    std::vector<unsigned long long> exp(bases.size(), 0);

                    for (size_t j = 0; j < bases.size(); j++)
                    {
                        if (possible_smooth[j].first == (i + offset) % bases[j] || possible_smooth[j].second == (i + offset) % bases[j])
                        {
                            while (smooth_residue % bases[j] == 0)
                            {
                                smooth_residue /= bases[j];
                                exp[j]++;
                            }
                        }
                    }

                    if (smooth_residue == 1)
                    {
                        #pragma omp critical
                        {
                            if (smooth_index.size() < bases.size() + 1)
                            {
                                smooth_index.push_back(smooth);
                                matrix.push_back(exp);

                                if (smooth_index.size() >= bases.size() + 1)
                                {
                                    done = true;
                                    #pragma omp flush(done)
                                }
                            }
                        }
                    }
                }
            }

            std::cout << "Thread " << omp_get_thread_num() << " terminou.\n";
        }

        offset += bases.size();
        std::cout << "Todas as threads terminaram.\n";
    }

    
    

    // #pragma omp barrier


    // omp_set_num_threads(1);


    // smooth_index.resize(bases.size() +1);
    // matrix.resize(bases.size() +1);
        

    
}

void quadratic_sieve(std::vector<mpz_class> &primes, mpz_class n) {

    std::vector<mpz_class> bases;
    bases.push_back(2);

    for (size_t i = 0; i < primes.size(); i++)
    {
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
        {
            bases.push_back(primes[i]);
        }
    }

    std::vector<std::pair<mpz_class, mpz_class>> smooth_index;
    smooth_index.reserve(bases.size() + 1);
    std::vector<std::vector<unsigned long long>> linear_system;
    linear_system.reserve(bases.size() + 1);

    sieving(bases, n, smooth_index, linear_system);
    
    // resolve the linear system                           


    std::vector<std::vector<int>> solutions = gauss_jordan(linear_system);


    mpz_class a, b;
    for(size_t combination=0; combination < solutions.size(); combination++) {       
        a = 1, b = 1;

        for (size_t i = 0; i < smooth_index.size(); i++) {
            if (solutions[combination][i]) {
                    
                a *= smooth_index[i].second;
                b *= smooth_index[i].first;
            }
        }

        a = std::sqrt(a);
        
        if((b - a) % n != 0 && (b + a) % n != 0)
        {
            std::cout << combination << '\n';
            break;
        }
    }

    std::cout << b << " " << a << '\n';
    std::cout << b - a << " " << a + b << '\n';
    std::cout << (b - a) % n << " " << (a + b) % n << '\n';
    mpz_class f1 = mdc(b - a, n);
    mpz_class f2 = mdc(a + b, n);
    std::cout << "n = " << f1 << " * " << f2 << "\n";

    
}