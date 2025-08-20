#pragma once
#include <vector>
#include <Eigen/Dense>

std::vector<Eigen::MatrixXd> applyMedianFilter(std::vector<Eigen::MatrixXd> channels, int windowSize);