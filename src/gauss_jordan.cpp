#include "gauss_jordan.hpp"

std::vector<int> gauss_jordan(std::vector<std::vector<unsigned long long>>& matrix) {
    int rows = matrix.size();
    int cols = matrix[0].size();

    for (int r = 0; r < rows; ++r) 
        for (int c = 0; c < cols; ++c)
            matrix[r][c] %= 2;

    std::vector<int> s(rows, 0); 
    for (int col = 0; col < cols; ++col) {
        int pivot_row = -1;
        for (int row = 0; row < rows; ++row) {
            if (matrix[row][col] == 1) {
                pivot_row = row;
                break;
            }
        }
        
        if (pivot_row == -1) continue; 
        
        s[pivot_row] = 1;

        for (int row = 0; row < rows; ++row) {
            if (row != pivot_row && matrix[row][col] == 1) {
                for (int c = 0; c < cols; ++c) {
                    matrix[row][c] ^= matrix[pivot_row][c]; 
                }
            }
        }
    }
    return s;
}