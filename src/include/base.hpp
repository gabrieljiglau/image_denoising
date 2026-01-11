#pragma once
#include <vector>
#include <Eigen/Dense>

class IAlgorithm{

    public: 
        
        virtual std::vector<std::vector<Eigen::MatrixXd>> apply();
        virtual ~IAlgorithm() = default;
};