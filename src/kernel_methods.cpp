#include "include/kernel_methods.hpp"
#include "include/utils.hpp"


Eigen::MatrixXd GaussianBlur::kernel(int kernelSize, double stddev){

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

Eigen::MatrixXd GaussianBlur::blur(Eigen::MatrixXd A, int kernelSize, double stddev){

    int padSize = kernelSize / 2;
    Eigen::MatrixXd kernel = GaussianBlur::kernel(kernelSize, stddev);
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

std::vector<Eigen::MatrixXd> GaussianBlur::applyBlur(int kernelSize, double stddev){
    
    std::vector<Eigen::MatrixXd> finalMatrices(3);
    std::vector<Eigen::MatrixXd> channels = getChannels();
    
    for (int i = 0; i < finalMatrices.size(); i++){
        // std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = blur(channels[i], kernelSize, stddev);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}

Eigen::MatrixXd MedianFilter::filter(Eigen::MatrixXd A, int windowSize){

    int padSize = windowSize / 2;

    Eigen::MatrixXd paddedMatrix = padMatrixReflect(A, padSize);
    Eigen::MatrixXd A_final = A;

    for (int row = 0; row < A.rows(); row++){
        for (int col = 0; col < A.cols(); col++){
            
            // std::cout << "Now at (row, col) = " << row << ", " << col << std::endl;
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

std::vector<Eigen::MatrixXd> MedianFilter::applyFilter(std::vector<Eigen::MatrixXd> channels, int windowSize){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        // std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = MedianFilter::filter(channels[i], windowSize);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}