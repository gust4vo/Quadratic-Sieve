#include "gauss_jordan.hpp"
#include <iostream>

std::vector<std::vector<int>> gauss_jordan(std::vector<std::vector<int>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    std::vector<std::vector<int>> gaus_reg(rows, std::vector<int>(rows, 0));

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

        
        if(pivot_row != pivot_qtd) {
            std::swap(matrix[pivot_qtd], matrix[pivot_row]);
            std::swap(indices[pivot_qtd], indices[pivot_row]);
        }

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

    std::vector<std::vector<int>> solution(rows - pivot_qtd);

    for (int i = 0; i < rows - pivot_qtd; i++)
        solution[i] = (gaus_reg[indices[i + pivot_qtd]]);

    return solution;
}