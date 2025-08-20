#include "medianFilter.hpp"
#include "utils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <algorithm>
#include <vector>
#include <iostream>
#include <Eigen/Dense>


static Eigen::MatrixXd medianFilter(Eigen::MatrixXd A, int windowSize){

    int padSize = windowSize / 2;

    Eigen::MatrixXd paddedMatrix = padMatrixReflect(A, padSize);
    Eigen::MatrixXd A_final = A;

    for (int row = 0; row < A.rows(); row++){
        for (int col = 0; col < A.cols(); col++){
            
            std::cout << "Now at (row, col) = " << row << ", " << col << std::endl;
            std::vector<double> window;

            int paddedRow = row + padSize;
            int paddedCol = col + padSize;

            // initialize window
            for (int i = -padSize; i <= padSize; i++){
                for (int j = -padSize; j <= padSize; j++){
                    window.push_back(paddedMatrix(paddedRow + i, paddedCol + j));
                }
            }

            std::sort(window.begin(), window.end());
            double medianValue = window[window.size() / 2];

            A_final(row,col) = medianValue;
        }
    }

    return A_final;
}

std::vector<Eigen::MatrixXd> applyMedianFilter(std::vector<Eigen::MatrixXd> channels, int windowSize){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = medianFilter(channels[i], windowSize);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}