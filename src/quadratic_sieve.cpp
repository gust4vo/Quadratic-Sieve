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
    
    std::vector<mpz_class> bases;
    bases.push_back(2);

    for (size_t i = 0; i < primes.size(); i++)
    {
        if (mpz_legendre(n.get_mpz_t(), primes[i].get_mpz_t()) == 1)
        {
            bases.push_back(primes[i]);
        }
    }

    mpz_class tonelli_ans, ans;
    size_t smooth_count = 0;
    std::vector<mpz_class> smooth_numbers(bases.size()+1), smooth_residue(bases.size()+1);
    mpz_class start = x; 
    
    while (smooth_count < (bases.size() + 1)) {
        std::vector<mpz_class> possible_smooth(primes.size()*6), sieve_list(primes.size()*6);

        for (size_t i = 0; i < primes.size()*6; i++) {
            possible_smooth[i] = (start + i); 
            sieve_list[i] = (possible_smooth[i] * possible_smooth[i]) % n;
        }

        tonelli_ans = cipolla(n, 2);
        ans = tonelli_ans - (start % 2);
        for (size_t k = ans.get_ui(); k < primes.size()*6; k += 2) {
            while (sieve_list[k] % 2 == 0) {
                sieve_list[k] /= 2;
            }
        }

        for (size_t i = 1; i < bases.size(); i++) {
            tonelli_ans = cipolla(n, bases[i]);
            ans = (tonelli_ans - (start % bases[i])) % bases[i];
            if (ans < 0) ans += bases[i];

            for (size_t k = ans.get_ui(); k < primes.size()*6; k += bases[i].get_ui()) {
                while (sieve_list[k] % bases[i] == 0) {
                    sieve_list[k] /= bases[i];
                }
            }

            ans = ((bases[i] - tonelli_ans) - (start % bases[i])) % bases[i];
            if (ans < 0) ans += bases[i];

            for (size_t k = ans.get_ui(); k < primes.size()*6; k += bases[i].get_ui()) {
                while (sieve_list[k] % bases[i] == 0) {
                    sieve_list[k] /= bases[i];
                }
            }
        }

        for (size_t i = 0; i < primes.size()*6; i++) {
            if (sieve_list[i] == 1) {
                smooth_numbers[smooth_count] = possible_smooth[i];
                smooth_residue[smooth_count] = (possible_smooth[i] * possible_smooth[i]) % n;
                smooth_count++;
                std::cout << smooth_count << " " << bases.size() + 1 << std::endl; 
                if (smooth_count == bases.size() + 1) break;
            }
        }

        start += primes.size()*6;
    }

    std::vector<std::vector<int>> exponents(smooth_count, std::vector<int>(bases.size(), 0));

    for(size_t i = 0; i < smooth_count; i++) {
        mpz_class aux = smooth_residue[i];
        for(size_t j = 0; j < bases.size(); j++) {
            while (aux % bases[j] == 0) {
                aux /= bases[j];
                exponents[i][j]++;
            }
        }
    }
    // resolve the linear system                           
    std::vector<std::vector<int>> solutions = gauss_jordan(exponents);
    mpz_class a = 1, b = 1;

    for(size_t combination=0; combination < solutions.size(); combination++) {
        for (size_t i = 0; i < smooth_count; i++) {
            if (solutions[combination][i]) {
                a *= smooth_residue[i];
                b *= smooth_numbers[i];
            }
        }

        mpz_sqrt(a.get_mpz_t(), a.get_mpz_t());

        if((b - a) % n != 0 && (a + b) % n != 0) break;

        a = 1, b = 1;
    }

    mpz_class f1 = mdc(b - a, n);
    mpz_class f2 = mdc(a + b, n);
    std::cout << "n = " << f1 << " * " << f2 << "\n";
}