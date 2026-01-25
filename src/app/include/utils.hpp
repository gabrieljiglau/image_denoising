#pragma once
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

int validateInput(std::string inputPng, bool checkExistence);

int addGaussianNoise(std::vector<unsigned char>image, int height, int width, int mean, std::vector<double> stdDevs, std::vector<std::string> savingPaths);

int addSaltPepperNoise(std::vector<unsigned char> image, std::string savingPath, int height, int width, float threshold);

Eigen::MatrixXd padMatrixReflect(Eigen::MatrixXd A, int padSize);

Eigen::VectorXd randomVector(int size, double epsilon);

std::vector<Eigen::MatrixXd> zeroMatrix(const int rows, std::vector<int> cols);

std::vector<Eigen::MatrixXd> zeroMatrix(const int rows, int col, int numDuplicates);

std::vector<Eigen::VectorXd> zeroVector(std::vector<int> cols);