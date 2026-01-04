#pragma once
#include <vector>
#include <Eigen/Dense>

std::vector<Eigen::MatrixXd> approximateImage(std::vector<Eigen::MatrixXd> channels, int k, int maxIterations, double epsilon);

std::vector<Eigen::MatrixXd> denoiseImage(std::vector<Eigen::MatrixXd> channels, int patchSize, int stride, int k, int maxIterations, double epsilon);

std::vector<Eigen::MatrixXd> applyMedianFilter(std::vector<Eigen::MatrixXd> channels, int windowSize);