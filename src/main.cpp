#include "lodepng.h"
#include "utils.hpp"
#include "powerIteration.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <vector>
#include <string>
#include <fmt/core.h>
#include <fmt/core.h>

int main(){

    std::vector<unsigned char> image;
    std::vector<unsigned char> newImage;
    unsigned width; // img_row
    unsigned height; // img_col

    std::string carPng = "/home/gabriel/Documents/HolyC/image_denoising/images/original.png";
    unsigned error = lodepng::decode(image, width, height, carPng);
    //printImage(image, height, width);

    //addNoise(image, height, width);

    std::vector<Eigen::MatrixXd> channels = rgbChannel(image, height, width);
    std::vector<int> kVals = {25, 50, 75, 100};

    int maxIterations = 10000;
    double epsilon = 1e-10;

    for (int k : kVals){
        std::vector<Eigen::MatrixXd> finalChannels = compressImage(channels, k, maxIterations, epsilon);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/compressed_k{}.png", k);
        reconstructImage(finalChannels, newImage, newPath, height, width);
    }

    // 15:08 start

    return 0;
}