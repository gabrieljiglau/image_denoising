#include "include/toolkit.hpp"
#include "include/base.hpp"
#include "include/kernel_methods.hpp"
#include "include/power_iteration.hpp"
#include "../libs/lodepng.h"
#include <Eigen/src/QR/ColPivHouseholderQR.h>
#include <iostream>
#include <memory>
#include <vector>
    

int Toolkit::loadPng(std::string pngPath){

    unsigned error = lodepng::decode(this->image, this->width, this->height, pngPath);
    return error;
}

/// remove from the input the k values that do not comply to: k < rank(A)
void Toolkit::checkRank(){

    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(this->channels[0]);
    int rank = qr.rank();

    int oldSize = kVals.size();

    kVals.erase(
        std::remove_if(kVals.begin(), kVals.end(), 
        [rank](int k) {return k > rank; }),
        kVals.end()
    );

}

/**
* @brief Deconstructs the image into 3(R/G/B) Eigen::MatrixXd matrices
*/
void Toolkit::initChannels(){

    loadPng(this->imagePath);

    std::vector<std::vector<int>> channelMatrices(3);

    for (unsigned int col = 0; col < height; col++){ // height        
        for (unsigned int row = 0; row < width; row++){ //width
            int idx = 4 * (col * width + row);
            channelMatrices[0].push_back((int)this->image[idx]);
            channelMatrices[1].push_back((int)this->image[idx + 1]);
            channelMatrices[2].push_back((int)this->image[idx + 2]);
        }
    }

    for (int i = 0; i < channelMatrices.size(); i++){
        std::vector<int> channel = channelMatrices[i];
        Eigen::Map<Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>> mappedMatrix(channel.data(), height, width);
        this->channels.push_back(mappedMatrix.cast<double>());
    }
    
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
    if (error != 0) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
        return error;
    } else {
        std::cout << "Processed image saved to " << newPath << std::endl;
    }

    return 0;
}

/**
 * @brief Applies the transformation based on the given 'mode' and saves the reconstructed image to 𝘯𝘦𝘸𝘗𝘢𝘵𝘩
 *
 */
std::vector<int> Toolkit::processPng(std::string inputPng, std::vector<std::string> newPaths){

    if (newPaths.size() != kVals.size() && !kVals.empty()){
        std::cout << "There must be 1:1 ration between kVals and the paths to store the image";
        return std::vector<int>();
    }

    std::unique_ptr<IAlgorithm> algorithm;
    std::vector<int> errCodes;
    std::vector<std::vector<Eigen::MatrixXd>> truncatedChannels;

    if (this->mode == "blur"){
        algorithm = std::make_unique<GaussianBlur>(channels, kernelSize, stddev);
    } else if (this->mode == "median_filter"){
        algorithm = std::make_unique<MedianFilter>(channels, windowSize);
    } else if (this->mode == "svd"){
        algorithm = std::make_unique<SVD>(channels, kVals, maxIterations, epsilon);
    } else if (this->mode == "patch_svd"){
        SVD patchSVD = SVD(channels, kVals ,maxIterations, epsilon);    
        truncatedChannels = patchSVD.apply(patchSize, stride);
    }

    if (algorithm){ // it wasn't patch svd
        truncatedChannels = algorithm->apply();
    }
    
    std::vector<std::vector<unsigned char>> newImages(newPaths.size());
    for (int i = 0; i < newPaths.size(); i++){
        errCodes.push_back(reconstructImage(truncatedChannels[i], newImages[i], newPaths[i]));
    }
    
    return errCodes;
}
