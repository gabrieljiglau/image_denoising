#include "lodepng.h"
#include "utils.hpp"
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

    std::string carPng = "/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k50.png";
    unsigned error = lodepng::decode(image, width, height, carPng);
    //printImage(image, height, width);
    
    //addNoise(image, height, width);
    std::vector<Eigen::MatrixXd> channels = rgbChannel(image, height, width);

    std::vector<int> kVals = {40};
    int maxIterations = 10000;
    double epsilon = 1e-10;

    // SVD after the median filter
    /*
    for (int k : kVals){
        std::vector<Eigen::MatrixXd> finalChannels = approximateImage(channels, k, maxIterations, epsilon);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k{}.png", k);
        reconstructImage(finalChannels, newImage, newPath, height, width);
    }
    */

    // de mutat fiecare flitru dinasta in propriul lui hpp

    /*
    // median filter (before/after SVD); still bad results
    int windowSize = 3;
    std::vector<Eigen::MatrixXd> medianChannels = applyMedianFilter(channels, windowSize);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_median_filter.png");
    reconstructImage(medianChannels, newImage, newPath, height, width);
    */

    // svd on smaller submatrices of the original image, bad results
    /*=
    int patchSize = 16;
    int stride = 12;
    int k = 2;
    
    std::vector<Eigen::MatrixXd> finalChannels = denoiseImage(channels, patchSize, stride, k, maxIterations, epsilon);
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
    
    // i) de reparat paddingu' atunci cand windows size > 3
    // ii) implementat gaussian blur

    return 0;
}