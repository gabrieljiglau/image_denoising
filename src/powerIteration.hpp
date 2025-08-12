#pragma once
#include <filesystem>
#include <fmt/core.h>
#include <vector>
#include <string>
#include <Eigen/Dense>

Eigen::MatrixXd powerIterationChannel(Eigen::MatrixXd A, int k, int maxIterations, double epsilon);

std::vector<Eigen::MatrixXd> compressImage(std::vector<Eigen::MatrixXd> channels, int k, int maxIterations, double epsilon);

