#pragma once
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

std::vector<Eigen::MatrixXd> rgbChannel(std::vector<unsigned char> image, int height, int width);

void addGaussianNoise(std::vector<unsigned char>image, int height, int width, std::vector<int> stdDevs);

void addSaltPepperNoise(std::vector<unsigned char> image, std::vector<unsigned char> &newImage, std::string newPath, 
    int height, int width, float threshold);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);

void reconstructImage(const std::vector<Eigen::MatrixXd> &truncatedChannels, std::vector<unsigned char> &newImage,
    std::string newPath, unsigned height, unsigned width);

void cumulateWeights(Eigen::MatrixXd &weights, int rowStart, int colStart, int patchSize);

Eigen::MatrixXd averagePixels(Eigen::MatrixXd A, Eigen::MatrixXd weights);

Eigen::MatrixXd padMatrixReflect(const Eigen::MatrixXd &A, int padSize);