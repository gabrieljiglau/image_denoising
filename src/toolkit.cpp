#include "include/toolkit.hpp"
#include "include/base.hpp"
#include "include/kernel_methods.hpp"
#include "libs/lodepng.h"
#include <iostream>


int Toolkit::loadPng(std::string pngPath){

    std::vector<unsigned char> newImage; // out image for reconstruction
    unsigned width; // img_row
    unsigned height; // img_col

    unsigned error = lodepng::decode(image, width, height, pngPath);
    return error;
}

/**
 * @brief Reconstructs an image (all 3 channels) and saves it to the given path
 *
 * @param truncatedChannels Each of the channels that constitute the image (R, G, B)
 * @param newImage Empty vector that will store the actual values of the pixels
 * @param newPath Where to save the image (PNG format only) 
 * @return The error message from encoding if that's the case, or 0 if the reconstruction went through
 *
 */
int Toolkit::reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage, std::string newPath){
    
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
        int a = 255; // I already know 'a' is 255 (the image isn't transparent)

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
        return error;
    }

    return 0;
}

int Toolkit::processPng(std::string inputPng){

    IAlgorithm algorithm;

    if (mode == "blur"){
        algorithm = GaussianBlur(channels, kernelSize, stddev);
    } else if (mode == "filter"){
        // median filter
    } else if (mode == "svd"){
        // svd
    } else if (mode == "patch_svd"){
        // patch svd
    }

    algorithm.apply();
    // reconstruct image;
    //save
}
