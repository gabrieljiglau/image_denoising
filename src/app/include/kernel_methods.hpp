#pragma once
#include "base.hpp"


class GaussianBlur: public IAlgorithm{

    protected:

        std::vector<Eigen::MatrixXd> channels;
        int kernelSize;
        double stddev;

        static Eigen::MatrixXd kernel(int kernelSize, double stddev);

        static Eigen::MatrixXd blur(Eigen::MatrixXd A, int kernelSize, double stddev);

    public:

        GaussianBlur(std::vector<Eigen::MatrixXd> channels, int kernelSize, double stddev): 
                    channels(channels), kernelSize(kernelSize), stddev(stddev){}

        std::vector<std::vector<Eigen::MatrixXd>> apply() override;

};

class MedianFilter: public IAlgorithm{

    protected:

        std::vector<Eigen::MatrixXd> channels;
        int windowSize;
        
        static Eigen::MatrixXd filter(Eigen::MatrixXd A, int windowSize);

    public:

        MedianFilter(std::vector<Eigen::MatrixXd> channels, int windowSize): channels(channels), windowSize(windowSize){}

        std::vector<std::vector<Eigen::MatrixXd>>  apply() override;

};