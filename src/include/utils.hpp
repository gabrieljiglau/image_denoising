#pragma once
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

void addGaussianNoise(std::vector<unsigned char>image, int height, int width, std::vector<int> stdDevs, std::string savingPath);

void addSaltPepperNoise(std::vector<unsigned char> image, std::vector<unsigned char> &newImage, std::string savingPath, 
                        int height, int width, float threshold);

void printImage(std::vector<unsigned char> imagePixels, int height, int width);

Eigen::MatrixXd padMatrixReflect(const Eigen::MatrixXd &A, int padSize);

Eigen::VectorXd randomVector(int size, double epsilon);