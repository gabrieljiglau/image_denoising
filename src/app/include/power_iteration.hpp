#pragma once
#include "base.hpp"
#include <vector>

class SVD: public IAlgorithm{

    protected:

        std::vector<Eigen::MatrixXd> channels;

        std::vector<int> kVals;
        int maxIterations;
        double epsilon;

        int patchSize;
        int stride;

        std::vector<Eigen::MatrixXd> powerIteration(Eigen::MatrixXd A, std::vector<int> k, int maxIterations, double epsilon);
        std::vector<Eigen::MatrixXd> applyPatches(Eigen::MatrixXd A, int patchSize, int stride, std::vector<int> k, int maxIterations, double epsilon);

        void cumulateWeights(Eigen::MatrixXd &weights, int rowStart, int colStart, int patchSize);
        std::vector<Eigen::MatrixXd> averagePixels(std::vector<Eigen::MatrixXd> A, std::vector<Eigen::MatrixXd> weights);


    public:

        SVD(std::vector<Eigen::MatrixXd> channels, std::vector<int> kVals, int maxIterations, double epsilon):
            channels(channels),
            kVals(kVals), 
            maxIterations(maxIterations), 
            epsilon(epsilon) {}

        virtual std::vector<std::vector<Eigen::MatrixXd>> apply() override;

        std::vector<std::vector<Eigen::MatrixXd>> apply(int patchSize, int stride);

};