#include "utils.hpp"
#include "lodepng.h"
#include <Eigen/src/Core/Matrix.h>
#include <random>
#include <cmath>
#include <iostream>
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

float gaussianGenerator(int mean, int stdDev){

    constexpr float PI = 3.14159265358979323846f;

    std::random_device rd;
    std::mt19937 generator(rd());

    // uniform distribution
    std::uniform_real_distribution<float> uniformDistribution(0.0f, 1.0f);

    float u1 = uniformDistribution(generator);
    float u2 = uniformDistribution(generator);

    // Box-Muller generator fod Gaussian pdf
    float z0 = std::sqrt(-2 * std::log(u1)) * std::cos(2 * PI * u2);

    return mean + stdDev * z0;
}

int addPixelNoise(int originalPixel, int mean, int stdDev){

    int noise = (int) gaussianGenerator(mean, stdDev);
    return (originalPixel + noise) % 256;

}

void generateNoisyImage(std::vector<unsigned char> originalImage, std::vector<unsigned char> &newImage, 
    std::string newPath, int height, int width, int mean, int stdDev){

    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){
            int idx = 4 * (col * width + row);

            int r = (int) originalImage[idx];
            int g = (int) originalImage[idx + 1];
            int b = (int) originalImage[idx + 2];
            int a = (int) originalImage[idx + 3];

            r += addPixelNoise(r, mean, stdDev);
            g += addPixelNoise(g, mean, stdDev);
            b += addPixelNoise(b, mean, stdDev);

            r = std::max(0, std::min(255, r));
            g = std::max(0, std::min(255, g));
            b = std::max(0, std::min(255, b));

            newImage.push_back(r);
            newImage.push_back(g);
            newImage.push_back(b);
            newImage.push_back(a);
        }
    }    

    int error = lodepng::encode(newPath, newImage, width, height);
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
    }

}

void addNoise(std::vector<unsigned char>image, int height, int width){

    int mean = 0; //mu
    std::vector<int> stdDevs = {2, 10};
    std::string newPath;
    std::vector<unsigned char> newImage;

    for (int stdDev : stdDevs){
        newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/noisy_mu{}_std{}.png", mean, stdDev);
        generateNoisyImage(image, newImage, newPath, height, width, mean, stdDev);
        std::cout << "Noisy image generated and saved to " << newPath << std::endl;
    }
}

std::vector<Eigen::MatrixXd> rgbChannel(std::vector<unsigned char> image, int height, int width){

    std::vector<std::vector<int>> channelMatrices(3);
    std::vector<Eigen::MatrixXd> returnVector(3);

    for (unsigned int col = 0; col < height; col++){ // height        
        for (unsigned int row = 0; row < width; row++){ //width
            int idx = 4 * (col * width + row);
            channelMatrices[0].push_back(image[idx]);
            channelMatrices[1].push_back(image[idx + 1]);
            channelMatrices[2].push_back(image[idx + 2]);
        }
    }

    for (int i = 0; i < returnVector.size(); i++){
        std::vector<int> channel = channelMatrices[i];
        Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channel.data(), height, width);
        returnVector[i] = mappedMatrix.cast<double>();
    }
    
    return returnVector;
}

void reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage,
    std::string newPath, unsigned height, unsigned width){
    
    std::vector<Eigen::RowVectorXd> rowVectors;

    for (Eigen::MatrixXd channel: truncatedChannels){
        Eigen::RowVectorXd innerVector(channel.size());

        // Eigen stores matrices in column-major order
        int idx = 0; 
        for (int row = 0; row < channel.rows(); row++){
            for (int col = 0; col < channel.cols(); col++){
                innerVector(idx++) = channel(row, col); 
            }
        }
        rowVectors.push_back(innerVector);
    }

    int size = rowVectors[0].size();
    for (int i = 0; i < size; i++){

        int r = (int) rowVectors[0](i);
        int g = (int) rowVectors[1](i);
        int b = (int) rowVectors[2](i);
        int a = 255; // I already know a is 255 (the image isn't transparent)

        r = std::max(0, std::min(255, r));
        g = std::max(0, std::min(255, g));
        b = std::max(0, std::min(255, b));

        newImage.push_back(r);
        newImage.push_back(g);
        newImage.push_back(b);
        newImage.push_back(a);
    }

    int error = lodepng::encode(newPath, newImage, width, height);
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
    }
}

Eigen::VectorXd randomVector(int size, double epsilon){
    
    Eigen::VectorXd vector = Eigen::VectorXd(size);

    std::normal_distribution<double> normalDistribution(0, 1);
    std::mt19937 generator(123456);

    for (int i = 0; i < size; i++){
        vector(i) = normalDistribution(generator);
    }

    double norm = vector.norm();
    if (norm < epsilon){
        vector = Eigen::VectorXd::Ones(size);
        return vector;
    }

    return vector / norm;
}

void cumulateWeights(Eigen::MatrixXd &weights, int rowStart, int colStart, int patchSize){

    for (int i = 0; i < patchSize; i++){
        for (int j = 0; j < patchSize; j++){
            weights(rowStart + i, colStart + j) += 1.0;
        }
    }
}

Eigen::MatrixXd averagePixels(Eigen::MatrixXd A, Eigen::MatrixXd weights){

    assert (A.rows() == weights.rows() && A.cols() == weights.cols());

    for (int i = 0; i < A.rows(); i++){
        for (int j = 0; j < A.cols(); j++){
            if (weights(i, j) > 0){
                A(i, j) /= weights(i, j);
            }
        }
    }

    return A;
}


void reflect(std::vector<std::vector<double>> &paddedVector, int padSize, bool left){

    int startIndex;
    int stopIndex;
    if (left){
        startIndex = 0;
        stopIndex = padSize;
    } else {
        startIndex = paddedVector[0].size() - padSize;
        stopIndex = paddedVector[0].size();
    }

    for (int i = 0; i < paddedVector.size(); i++){
        std::vector<double> &currentVector = paddedVector[i];

        int reflectionIdx = currentVector.size() - 2 * padSize;
        for (int idx = startIndex; idx < stopIndex; idx++){
            currentVector[idx] = currentVector[reflectionIdx];
            reflectionIdx--;
        }
    }

}

std::vector<std::vector<double>> buildRowPad(Eigen::MatrixXd A, int padSize, bool top){

    std::vector<std::vector<double>> paddedVector;
    int startRow;
    if (top){
        startRow = padSize;
    } else {
        startRow = A.rows() - padSize - 1; // era -2
    }

    for (int padIdx = 0; padIdx < padSize; padIdx++){
        std::vector<double> innerPad(A.row(0).size() + 2 * padSize);
        for (int i = 0; i < A.row(0).size(); i++){
            innerPad[padSize + i] = A.row(startRow)(i);
        }

        startRow--;
        paddedVector.push_back(innerPad);
    }

    bool left = true;
    reflect(paddedVector, padSize, left);
    reflect(paddedVector, padSize, !left);
    return paddedVector;
}

std::vector<std::vector<double>> buildColumnPad(Eigen::MatrixXd A, int padSize, bool left){

    std::vector<std::vector<double>> paddedVector;
    int startCol;
    if (left){
        startCol = padSize;
    } else if (!left){
        startCol = A.cols() - padSize - 1; // era -2
    }

    for (int padIdx = 0; padIdx < padSize; padIdx++){

        std::vector<double> innerPad(A.col(0).size());
        for (int j = 0; j < A.col(0).size(); j++){
            innerPad[j] = A.col(startCol)(j);
        }

        startCol--;
        paddedVector.push_back(innerPad);
    }

    return paddedVector;
}

Eigen::MatrixXd padMatrix(Eigen::MatrixXd A, int padSize){
    
    // padSize * 2, since the pad is above the first row, and below the last row (and similar for the columns)
    Eigen::MatrixXd newA = Eigen::MatrixXd::Zero(A.rows() + padSize * 2, A.cols() + padSize * 2);
    newA.block(padSize, padSize, A.rows(), A.cols()) = A;

    // fill top padding
    std::vector<std::vector<double>> topPad = buildRowPad(A, padSize, /*top=*/true);
    for (int row = 0; row < topPad.size(); row++){
        for (int col = 0; col < newA.cols(); col++){
            newA(row, col) = topPad[row][col];
        }
    }

    // fill bottom padding
    std::vector<std::vector<double>> bottomPad = buildRowPad(A, padSize, /*top=*/false);
    for (int row = 0; row < padSize; row++){
        for (int col = 0; col < newA.cols();col++){
            newA(newA.rows() - padSize + row, col) = bottomPad[row][col];
        }
    }
    
    std::vector<std::vector<double>> leftPad = buildColumnPad(A, padSize, /*left=*/true);

    for (int i = 0; i < padSize; i++){
        std::vector<double> &innerPad = leftPad[i];
        for (int j = 0; j < innerPad.size(); j++){
            double element = innerPad[j];
            newA(padSize + j, i) = innerPad[j];
        }
    }  

    std::vector<std::vector<double>> rightPad = buildColumnPad(A, padSize, /*left=*/false);
    for (int i = 0; i < padSize; i++){
        std::vector<double> &innerPad = rightPad[i];
        for (int j = 0; j < innerPad.size(); j++){
            double element = innerPad[j];
            newA(padSize + j, newA.cols() - padSize + i) = innerPad[j];
        }
    }

    return newA;
}

