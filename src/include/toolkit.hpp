#pragma once
#include "base.hpp"
#include <string>
#include <vector>
#include <Eigen/Dense>

class Toolkit{

    private:

        std::vector<unsigned char> image;
        int width;
        int height;

        std::string mode;

        std::vector<Eigen::MatrixXd> channels;

        // for Gaussian Blur
        int kernelSize;
        double stddev;

        // for median filter
        int windowSize;

        // for svd
        std::vector<int> kVals;
        int maxIterations;
        double epsilon;

        // extra parametes for patch_svd
        int patchSize;
        int stride;

        IAlgorithm *algorithm;

        std::vector<Eigen::MatrixXd> rgbChannel();

    public:

        Toolkit(std::string imagePath, std::string mode, int width = 0, int height = 0): mode(mode), width(width), height(height) {    
            loadPng(imagePath);     // fills internally this->image
            this->channels = rgbChannel();  // transform the image into Eigen::MatrixXd 
        };
                

        int loadPng(std::string pngPath);

        int reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage, std::string newPath);

        //std::vector<Eigen::MatrixXd> getChannels() {return this->channels;}

        void checkRank();

        void setKernelSize(int kernelSize) {this->kernelSize = kernelSize;}
        void setStddev(double stddev) {this->stddev = stddev;}

        void setWindowSize(int windowSize) {this->windowSize = windowSize;}

        void setK(std::vector<int> kVals) {this->kVals = kVals;}
        void setMaxIterations(int maxIterations) {this->maxIterations = maxIterations;}
        void setEpsilon(double epsilon) {this->epsilon = epsilon;}

        void setPatchSize(int patchSize) {this->patchSize = patchSize;}
        void setStride(int stride) {this->stride = stride;}

        std::vector<int> processPng(std::string inputPng, std::string newPath);

};