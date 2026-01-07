#include "include/kernel_methods.hpp"
#include "include/utils.hpp"

/**
 * @brief Prepares the gaussian kernel to be applied on a given channel
 *
 * @param kernelSize  The size of the kernel (usually an odd number)
 * @param stddev The standard deviation used in gaussian distribution ~ N (0, stddev).
 * @return The kernel as an Eigen::MatrixXd.
 *
 */
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

/// @param A: the channel for which the transformation is applied
/// @return sum of the multiplications between the original pixel and the elements from the kernel
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

/// Applies the Gaussian Blur on all 3 channels from the image
std::vector<Eigen::MatrixXd> GaussianBlur::apply(){
    
    std::vector<Eigen::MatrixXd> finalMatrices(3);
    
    for (int i = 0; i < finalMatrices.size(); i++){
        // std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = blur(channels[i], kernelSize, stddev);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}


/**
 * @brief Prepares the median filter to be applied on a given channel
 *
 * @param A The channel for which the transformation is applied
 * @param windowSize  The size of the kernel (usually an odd number)
 * @return The filter as an Eigen::MatrixXd.
 *
 */
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


/// Applies the Median Filter on all 3 channels from the image
std::vector<Eigen::MatrixXd> MedianFilter::apply(){

    std::vector<Eigen::MatrixXd> finalMatrices(3);

    for (int i = 0; i < finalMatrices.size(); i++){
        // std::cout << "Now at index = " << i << std::endl;
        Eigen::MatrixXd newChannel = MedianFilter::filter(channels[i], windowSize);
        finalMatrices[i] = newChannel;
    }

    return finalMatrices;
}