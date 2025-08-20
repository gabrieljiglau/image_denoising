#include "lodepng.h"
#include "utils.hpp"
#include "medianFilter.hpp"
#include "gaussianBlur.hpp"
#include "powerIteration.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/core.h>
#include <fmt/core.h>

int main(){

    std::vector<unsigned char> image;

    // out for reconstruction
    std::vector<unsigned char> newImage;
    unsigned width; // img_row
    unsigned height; // img_col

    
    std::string carPng = "/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k50.png";
    unsigned error = lodepng::decode(image, width, height, carPng);
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
    /*
    for (int k : kVals){
        std::vector<Eigen::MatrixXd> finalChannels = approximateImage(channels, k, maxIterations, epsilon);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k{}.png", k);
        reconstructImage(finalChannels, newImage, newPath, height, width);
    }
    */

    // gaussianBlur after SVD
    int kernelSize = 3;
    std::vector<int> stdDevs = {2, 5, 10, 20};
    for (int stdDev : stdDevs){
        std::vector<Eigen::MatrixXd> blurredChannels = applyGaussianBlur(channels, kernelSize, stdDev);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k{}_afterBlur.png", stdDev);
        reconstructImage(blurredChannels, newImage, newPath, height, width);
    }



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