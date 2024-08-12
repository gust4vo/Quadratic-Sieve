#include "quadratic_sieve.hpp"
#include "tonelli_shanks.hpp"
#include "cipolla.hpp"
#include "mdc.hpp"
#include "ak_mod_n.hpp"
#include "gauss_jordan.hpp"
#include <cmath>
#include <iostream>

void quadratic_sieve(std::vector<mpz_class> &primes, mpz_class n) {

    mpz_class x = (std::sqrt(n));

    if (x*x < n) x++;

    size_t possible_smooth_size = primes.size() * 6;

    std::vector<mpz_class> possible_smooth(possible_smooth_size);

    for (size_t i = 0; i < possible_smooth_size; i++) {
        possible_smooth[i] = (x + i) * (x + i) % n;
    }
    
    std::vector<mpz_class> bases;
    bases.push_back(2);

    std::cout << "\nBase: ";

    for (size_t i = 0; i < primes.size(); i++)
    {
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
        {
            bases.push_back(primes[i]);
            std::cout << primes[i] << ' ';
        }
    }
        std::cout << '\n';

    std::vector<std::vector<unsigned long long>> pre_matriz(
    possible_smooth_size,
    std::vector<unsigned long long>(bases.size(), 0));

    mpz_class tonelli_ans, ans;
    
    tonelli_ans = cipolla(n, 2);

    ans = tonelli_ans - (x % 2);
    for (size_t k = ans.get_ui(); k < possible_smooth_size; k += 2) {
        // std::cout << k << " " << possible_smooth[k] << std::endl; 
        while(possible_smooth[k] % 2 == 0) 
        {
            possible_smooth[k] /= 2;
            pre_matriz[k][0]++;
        }
    }

    // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[0] << '\n';

    // for (size_t k = 0; k < possible_smooth_size; k++)
    //     std::cout << possible_smooth[k] << ' ';

    // std::cout << '\n';
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

        // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[i] << '\n';
        
        // for (size_t k = 0; k < possible_smooth_size; k++)
        //     std::cout << possible_smooth[k] << ' ';

        // std::cout << '\n';
    
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

        // std::cout << "Toneli: " << tonelli_ans << " Begin: " << ans << " Prime :" << bases[i] << '\n';

        // for (size_t k = 0; k < possible_smooth_size; k++)
        //     std::cout << possible_smooth[k] << ' ';

        // std::cout << '\n';
        
        
    }

    // for (size_t i = 0; i < possible_smooth_size; i++)
    // {
    //     std::cout << '\t' << (x + i)*(x + i) << "\t"  << (x + i)*(x + i) % n << "\t" << possible_smooth[i] << "\t\t";

    //     for (size_t j = 0; j < bases.size(); j++)
    //     {
    //         std::cout << pre_matriz[i][j] << ' ';
    //     }

    //     std::cout << '\n';
    // }

    std::vector<size_t> smooth_index;
    for (size_t i = 0, aux = 0; i < possible_smooth_size; i++) {
        if (possible_smooth[i] == 1) {
            smooth_index.push_back(i);
            aux++;
        }

        if (aux == bases.size() + 1) break;
    }

    if (smooth_index.size() < bases.size() + 1) {
        std::cout << "There's no sufficient smooth numbers" << std::endl;
        return;
    }

    std::vector<std::vector<unsigned long long>> linear_system;
    for (auto index : smooth_index) {
        linear_system.push_back(pre_matriz[index]);
    }
    
    // resolve the linear system                           
    std::vector<std::vector<int>> solutions = gauss_jordan(linear_system);

    std::cout <<  '\n';

    mpz_class a = 1, b = 1, aux = 1;
    for(unsigned long int combinacao=0; combinacao < solutions.size(); combinacao++) {
        std::cout << "Combinacao " << combinacao << " -> ";
        for (size_t i = 0; i < smooth_index.size(); i++) {
            if (solutions[combinacao][i]) {
                for(size_t j = 0; j < bases.size(); j++) {
                    mpz_class p;

                    // std::cout << pre_matriz[smooth_index[i]][j] << " ";

                    mpz_pow_ui(p.get_mpz_t(), bases[j].get_mpz_t(), pre_matriz[smooth_index[i]][j]);
                    a *= p;
                    aux *= p;
                }

                // std::cout << " -> " << aux << "\n";
                aux = 1;
                b *= x + smooth_index[i];
            }
        }

        mpz_sqrt(a.get_mpz_t(), a.get_mpz_t());

        // std::cout << "a = " << a << " b = " << b  << '\n';
        mpz_class f1 = mdc(b - a, n);
        mpz_class f2 = mdc(a + b, n);
        std::cout << "n = " << f1 << " * " << f2 << "\n";
        a = 1, b = 1;
    }
}