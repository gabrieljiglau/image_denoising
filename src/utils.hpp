#pragma once
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

std::vector<Eigen::MatrixXd> rgbChannel(std::vector<unsigned char> image, int height, int width);

void addNoise(std::vector<unsigned char>image, int height, int width);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);

void reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage,
    std::string newPath, unsigned height, unsigned width);