#include "include/toolkit.hpp"
#include "include/base.hpp"
#include "include/kernel_methods.hpp"
#include "include/power_iteration.hpp"
#include "libs/lodepng.h"
#include <Eigen/src/QR/ColPivHouseholderQR.h>
#include <iostream>
    

int Toolkit::loadPng(std::string pngPath){

    unsigned error = lodepng::decode(this->image, this->width, this->height, pngPath);
    return error;
}

/// remove from the input the k values that do not comply to: k < rank(A)
void Toolkit::checkRank(){

    /// TODO: functia asta 'sparta' intr-o functie de initializare si in una de verificare a k-ului 
    
    rgbChannel();
    std::cout << "channels[0].size()" << this->channels[1].size() << std::endl;
    Eigen::ColPivHouseholderQR<Eigen::MatrixXd> qr(this->channels[0]);
    int rank = qr.rank();

    int oldSize = kVals.size();

    kVals.erase(
        std::remove_if(kVals.begin(), kVals.end(), 
        [rank](int k) {return k > rank; }),
        kVals.end()
    );

}


void Toolkit::rgbChannel(){

    loadPng(this->imagePath);

    std::cout << "this.width = " << this->width << " this.height = " << this->height << std::endl;

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
    if (error) {
        std::cerr << "Encoder error: " << error << lodepng_error_text(error) << std::endl;
        return error;
    }

    return 0;
}

/**
 * @brief Applies the transformation based on the given 'mode' and saves the reconstructed image to 𝘯𝘦𝘸𝘗𝘢𝘵𝘩
 *
 */
std::vector<int> Toolkit::processPng(std::string inputPng, std::vector<std::string> newPaths){

    IAlgorithm *algorithm;
    std::vector<int> errCodes;
    std::vector<std::vector<Eigen::MatrixXd>> truncatedChannels;

    if (this->mode == "blur"){
        algorithm = new GaussianBlur(channels, kernelSize, stddev);
    } else if (this->mode == "filter"){
        algorithm = new MedianFilter(channels, windowSize);
    } else if (this->mode == "svd"){
        algorithm = new SVD(channels, kVals, maxIterations, epsilon);
    } else if (this->mode == "patch_svd"){
        SVD patchSVD = SVD(channels, kVals ,maxIterations, epsilon);
        truncatedChannels = patchSVD.apply(patchSize, stride);
    }

    if (truncatedChannels.empty()){ // it wasn't patch svd
        truncatedChannels = algorithm->apply();
    }
    
    std::vector<unsigned char> newImage;
    for (int i = 0; i < truncatedChannels.size(); i++){
        std::vector<Eigen::MatrixXd> channel = truncatedChannels[i];
        std::string newPath = newPaths[i];
        errCodes.push_back(reconstructImage(channel, newImage, newPath));
    }
    
    delete algorithm;
    algorithm = nullptr;
    
    return errCodes;
}
