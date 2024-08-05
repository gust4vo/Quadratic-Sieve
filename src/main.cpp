#include <iostream>
#include "tonelli_shanks.hpp"


int main()
{
    mpz_class n, p;
    std::cin >> n >> p;

    std::cout << tonelli_shanks(n, p) << '\n';

}