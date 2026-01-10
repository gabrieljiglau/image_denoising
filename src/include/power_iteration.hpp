#pragma once
#include "base.hpp"
#include <vector>

class SVD: public IAlgorithm{

    protected:

        std::vector<Eigen::MatrixXd> channels;

        int k;
        int maxIterations;
        double epsilon;

        int patchSize;
        int stride;

        Eigen::MatrixXd powerIteration(Eigen::MatrixXd A, int k, int maxIterations, double epsilon);
        Eigen::MatrixXd applyPatches(Eigen::MatrixXd A, int patchSize, int stride, int k, int maxIterations, double epsilon);

        void cumulateWeights(Eigen::MatrixXd &weights, int rowStart, int colStart, int patchSize);
        Eigen::MatrixXd  averagePixels(Eigen::MatrixXd A, Eigen::MatrixXd weights);


    public:

        SVD(std::vector<Eigen::MatrixXd> channels, int k, int maxIterations, double epsilon):
            channels(channels),
            k(k), 
            maxIterations(maxIterations), 
            epsilon(epsilon) {}

        virtual std::vector<Eigen::MatrixXd> apply() override;

        std::vector<Eigen::MatrixXd> apply(int patchSize, int stride);

};