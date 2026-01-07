#include "include/lodepng.h"
#include "include/utils.hpp"
#include "include/medianFilter.hpp"
#include "include/gaussianBlur.hpp"
#include "include/toolkit.hpp"
#include "include/powerIteration.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <vector>
#include <string>
#include <optional>
#include <stdio.h>
#include <fmt/core.h>
#include <fmt/core.h>


int main(int argc, char *argv[]){

    if (argc != 3 || argc != 4) {
        printf("Usage: ./build/image_denoising full/path/to/image mode kForSVD \n");
        return 1;
    }

    std::string inputPng = argv[1];
    std::string mode = argv[2];
    std::optional<std::string> kSVD = argv[3];

    Toolkit toolkit;
    toolkit.processPng(inputPng, mode, kSVD);

    //printImage(image, height, width);

    /*
    std::vector<int> stdDevs = {10};
    addGaussianNoise(image, height, width, stdDevs);
    */
    /*
    addSaltPepperNoise(image, newImage, "/home/gabriel/Documents/HolyC/image_denoising/images/noisy/car_salt_pepper.png",
                     height, width, 0.1);
    */
    /*
    std::vector<Eigen::MatrixXd> channels = rgbChannel(image, height, width);
    int windowSize = 3;
    std::vector<Eigen::MatrixXd> medianChannels = applyMedianFilter(channels, windowSize);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_car_median.png");
    reconstructImage(medianChannels, newImage, newPath, height, width);
    */  
    
    
    std::vector<Eigen::MatrixXd> channels = rgbChannel(image, height, width);

    std::vector<int> kVals = {50};
    int maxIterations = 10000;
    double epsilon = 1e-10;


    // SVD
    for (int k : kVals){
        std::vector<Eigen::MatrixXd> finalChannels = approximateImage(channels, k, maxIterations, epsilon);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k{}.png", k);
        reconstructImage(finalChannels, newImage, newPath, height, width);
    }
    
    // gaussianBlur after SVD
    /*
    int kernelSize = 3;
    std::vector<int> stdDevs = {2, 5, 10, 20};
    for (int stdDev : stdDevs){
        std::vector<Eigen::MatrixXd> blurredChannels = applyGaussianBlur(channels, kernelSize, stdDev);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k{}_afterBlur.png", stdDev);
        reconstructImage(blurredChannels, newImage, newPath, height, width);
    }
    */



    /*
    // median filter (before/after SVD); still bad results
    int windowSize = 2;
    std::vector<Eigen::MatrixXd> medianChannels = applyMedianFilter(channels, windowSize);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k50_median.png");
    reconstructImage(medianChannels, newImage, newPath, height, width);
    */

    // svd on smaller submatrices of the original image, bad results
    /*=
    int patchSize = 16;
    int stride = 12;
    int k = 2;
    
    std::vector<Eigen::MatrixXd> finalChannels = patchesDecomp(channels, patchSize, stride, k, maxIterations, epsilon);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k{}.png", k);
    reconstructImage(finalChannels, newImage, newPath, height, width);
    */


    // compression (SVD on the original image)
    /*
    for (int k : kVals){
        std::vector<Eigen::MatrixXd> finalChannels = approximateImage(channels, k, maxIterations, epsilon);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/compressed_k{}.png", k);
        reconstructImage(finalChannels, newImage, newPath, height, width);
    }
    */

    return 0;
}