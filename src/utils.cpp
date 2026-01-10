#include "include/utils.hpp"
#include "libs/lodepng.h"
#include <Eigen/src/Core/Matrix.h>
#include <random>
#include <cmath>
#include <iostream>
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>


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

void saltPepperNoise(int &r, int &g, int &b, float threshold){

    std::random_device rd;
    std::mt19937 generator(rd());

    std::uniform_real_distribution<float> uniformDistribution(0.0f, 1.0f);

    double generatedValue = uniformDistribution(generator);
    double saltOrPepper = uniformDistribution(generator);
    if (generatedValue < threshold){
        if (saltOrPepper < 0.5){ // salt
            r = g = b = 0;
        } else { //pepper
            r = g = b = 255; 
        }
    }
}

int addPixelNoise(int originalPixel, int mean, int stdDev){

    int noise = (int) gaussianGenerator(mean, stdDev);
    return (originalPixel + noise) % 255;

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

void addSaltPepperNoise(std::vector<unsigned char> image, std::vector<unsigned char> &newImage, std::string newPath, 
                        int height, int width, float threshold){

    for (unsigned int col = 0; col < height; col++){
        for (unsigned int row = 0; row < width; row++){

            int idx = 4 * (col * width + row);

            int r = (int) image[idx];
            int g = (int) image[idx + 1];
            int b = (int) image[idx + 2];
            int a = (int) image[idx + 3];

            saltPepperNoise(r, g, b, threshold);

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

    std::cout << "Image with salt and pepper noise generated and saved to " << newPath << std::endl;
}


/// TODO: functia asta sa primeasca ca parametru locul de salvare al imaginii
void addGaussianNoise(std::vector<unsigned char>image, int height, int width, std::vector<int> stdDevs){

    int mean = 0; //mu
    std::string newPath;
    std::vector<unsigned char> newImage;

    for (int stdDev : stdDevs){
        newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/noisy/noisy_mu{}_std{}.png", mean, stdDev);
        generateNoisyImage(image, newImage, newPath, height, width, mean, stdDev);
        std::cout << "Image with noise (from a gaussian PDF) generated and saved to " << newPath << std::endl;
    }
}

Eigen::MatrixXd padMatrixReflect(const Eigen::MatrixXd &A, int padSize) {

    int rows = A.rows();
    int cols = A.cols();

    int newRows = rows + 2 * padSize;
    int newCols = cols + 2 * padSize;

    Eigen::MatrixXd A_final(newRows, newCols);

    // center
    A_final.block(padSize, padSize, rows, cols) = A;

    for (int i = 0; i < padSize; i++) {
        A_final.row(padSize - 1 - i) = A_final.row(padSize + i);              // reflect top
        A_final.row(padSize + rows + i) = A_final.row(padSize + rows - 1 - i); // reflect bottom
    }

    for (int i = 0; i < padSize; i++) {
        A_final.col(padSize - 1 - i) = A_final.col(padSize + i);               // reflect left
        A_final.col(padSize + cols + i) = A_final.col(padSize + cols - 1 - i); // reflect right
    }

    return A_final;
}
