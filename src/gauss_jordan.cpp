#include "gauss_jordan.hpp"
#include <iostream>

std::vector<std::vector<int>> gauss_jordan(std::vector<std::vector<unsigned long long>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int r = 0; r < rows; r++) 
        for (int c = 0; c < cols; c++)
            matrix[r][c] %= 2;

    std::vector<std::vector<int>> gaus_reg(rows, std::vector<int>(rows, 0));
    std::vector<std::vector<int>> solution;

    int *indices = new int[rows];

    for (int i = 0; i < rows; i++)
    {
        indices[i] = i;
        gaus_reg[i][i] = 1;
    }

    int pivot_qtd = 0;

    for (int col = 0; col < cols; col++) {
        int pivot_row = -1;
        for (int row = pivot_qtd; row < rows; row++) {
            if (matrix[row][col] == 1) {
                pivot_row = row;
                break;
            }
        }
        
        if (pivot_row == -1)
            continue; 

        // for (size_t i = 0; i < rows; i++)
        // {
        //     std::cout << indices[i] << "\t";
        //     for (size_t j = 0; j < cols; j++)
        //     {
        //         std::cout << matrix[i][j] << " ";
        //     }

        //     std::cout << '\n';
            
        // }

        // std::cout << '\n';
        
        if(pivot_row != pivot_qtd) {
            std::swap(matrix[pivot_qtd], matrix[pivot_row]);
            std::swap(indices[pivot_qtd], indices[pivot_row]);
        }

        // for (size_t i = 0; i < rows; i++)
        // {
        //     std::cout << indices[i] << "\t";
        //     for (size_t j = 0; j < cols; j++)
        //     {
        //         std::cout << matrix[i][j] << " ";
        //     }

        //     std::cout << '\n';
            
        // }

        // std::cout << '\n';

        for (int row = pivot_qtd + 1; row < rows; row++) {
            if (matrix[row][col] == 1) {                

                for (int c = 0; c < cols; c++) {
                    matrix[row][c] ^= matrix[pivot_qtd][c]; 
                    gaus_reg[indices[row]][c] ^= gaus_reg[indices[pivot_qtd]][c];

                }
            }
        }

        pivot_qtd++;
    }


    for (int i = 0; i < rows; i++) {
        bool is_null_row = true;
        for (int j = 0; j < cols; j++) {
            if (matrix[i][j] != 0) {
                is_null_row = false;
                break;
            }
        }
        if (is_null_row) {
            solution.push_back(gaus_reg[indices[i]]);
        }  
    }
    
    return solution;
}