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
        int k;
        int maxIterations;
        double epsilon;

        int patchSize;
        int stride;

        IAlgorithm algorithm;

        std::vector<Eigen::MatrixXd> rgbChannel();

    public:

        Toolkit(std::vector<unsigned char> image, std::string mode, int k, int width = 0, int height = 0):
            image(image), mode(mode), k(k), width(width), height(height) {};

        int loadPng(std::string pngPath);

        int processPng(std::string inputPng, std::string newPath);

        int reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage, std::string newPath);

        //std::vector<Eigen::MatrixXd> getChannels() {return this->channels;}

        void setKernelSize(int kernelSize) {this->kernelSize = kernelSize;}
        void setStddev(double stddev) {this->stddev = stddev;}

        void setWindowSize(int windowSize) {this->windowSize = windowSize;}

        void setK(int k) {this->k = k;}

        void setPatchSize(int patchSize) {this->patchSize = patchSize;}
        void setStride(int stride) {this->stride = stride;}

};