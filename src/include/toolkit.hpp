#pragma once
#include "base.hpp"
#include <optional>
#include <string>
#include <vector>
#include <Eigen/Dense>

class Toolkit{

    private:

        std::vector<unsigned char> image;
        int width;
        int height;

        std::string mode;
        std::optional<std::string> kSVD;

        std::vector<Eigen::MatrixXd> channels;

        // for Gaussian Blur
        int kernelSize;
        int stddev;

        IAlgorithm algorithm;

    public:

        Toolkit(std::vector<unsigned char> image, std::string mode, std::optional<std::string> kSVD, int witdth = 0, int height = 0):
            image(image),
            mode(mode),
            kSVD(kSVD) {};

        //std::vector<Eigen::MatrixXd> getChannels() {return this->channels;}

        int loadPng(std::string pngPath);

        int processPng(std::string inputPng);

        int reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage, std::string newPath);

};