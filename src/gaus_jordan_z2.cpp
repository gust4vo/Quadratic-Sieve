#include "gauss_jordan_z2.hpp"
#include <vector>
#include <iostream>
#include <stdexcept>

void print_matrix(const std::vector<std::vector<int>>& mat) {
    for (const auto& row : mat) {
        for (int val : row) {
            std::cout << val << " ";
        }
        std::cout << "\n";
    }
}

std::vector<int> gauss_jordan_z2(const std::vector<std::vector<int>>& matrix, const std::vector<int>& b) {
    int n = matrix.size();
    int m = matrix[0].size();

    if (n != b.size()) {
        throw std::invalid_argument("A dimensão da matriz não corresponde ao vetor b.");
    }

    std::vector<std::vector<int>> augmented(n, std::vector<int>(m + 1, 0));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            augmented[i][j] = matrix[i][j];
        }
        augmented[i][m] = b[i];
    }

    for (int col = 0; col < m; ++col) {
        int row = -1;
        for (int i = col; i < n; ++i) {
            if (augmented[i][col] == 1) {
                row = i;
                break;
            }
        }

        if (row == -1) continue;

        if (row != col) {
            std::swap(augmented[row], augmented[col]);
        }

        for (int i = col + 1; i < n; ++i) {
            if (augmented[i][col] == 1) {
                for (int j = col; j <= m; ++j) {
                    augmented[i][j] ^= augmented[col][j];
                }
            }
        }
    }

    std::vector<int> solution(m, 0);
    for (int i = m - 1; i >= 0; --i) {
        if (augmented[i][i] == 1) {
            solution[i] = augmented[i][m];
            for (int j = i + 1; j < m; ++j) {
                solution[i] ^= (augmented[i][j] * solution[j]);
            }
        }
    }

    return solution;
}
