#include "gaussianBlur.hpp"
#include "utils.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <vector>
#include <iostream>
#include <Eigen/Dense>

Eigen::MatrixXd gaussianKernel(int kernelSize, double stddev){

    Eigen::MatrixXd kernel(kernelSize, kernelSize);
    int half = kernelSize / 2;
    double sum = 0;

    // around the center
    for (int i = -half; i <= half; i++){
        for (int j = -half; j <= half; j++){
            double value = std::exp((i * i + j * j) / (-2 * stddev * stddev)); 
            kernel(i + half, j + half) = value;
            sum += value;
        }
    }

    return kernel / sum;
}

static Eigen::MatrixXd gaussianBlur(Eigen::MatrixXd A, int kernelSize, double stddev){

    int padSize = kernelSize / 2;
    Eigen::MatrixXd kernel = gaussianKernel(kernelSize, stddev);
    Eigen::MatrixXd paddedMatrix = padMatrixReflect(A, padSize);

    Eigen::MatrixXd A_final = A;
    for (int row = 0; row < A.rows(); row++){
        for (int col = 0; col < A.cols(); col++){
            
            double sum = 0;
            for (int i = 0; i < kernelSize; i++){
                for (int j = 0; j < kernelSize; j++){
                    sum += paddedMatrix(row + i, col + j) * kernel(i, j);
                }
            }

            A_final(row, col) = sum;
        }
    }

    return A_final;
}

std::vector<Eigen::MatrixXd> applyGaussianBlur(std::vector<Eigen::MatrixXd> channels, int kernelSize, double stddev){
    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = gaussianBlur(channels[i], kernelSize, stddev);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}