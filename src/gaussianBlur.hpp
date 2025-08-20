#pragma once
#include <vector>
#include <Eigen/Dense>

Eigen::MatrixXd gaussianKernel(int kernelSize, double stddev);

std::vector<Eigen::MatrixXd> applyGaussianBlur(std::vector<Eigen::MatrixXd> channels, int kernelSize, double stddev);