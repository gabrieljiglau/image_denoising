#pragma once
#include <vector>
#include <Eigen/Dense>

class IKernel{

    /// TODO: ai interfata asta, care va fi redenumita IAlgorithm, care va fi implementata de toate functiile aplicabile

    protected:
        std::vector<Eigen::MatrixXd> channels;

    public: 

        virtual std::vector<Eigen::MatrixXd> apply();
        virtual ~IKernel() = default;
};

class GaussianBlur: public IKernel{

    private:

        int kernelSize;
        double stddev;
        static Eigen::MatrixXd kernel(int kernelSize, double stddev);
        static Eigen::MatrixXd blur(Eigen::MatrixXd A, int kernelSize, double stddev);

    public:

        GaussianBlur(std::vector<Eigen::MatrixXd> channels, int kernelSize, double stddev) {
            this->channels = channels;
            this->kernelSize = kernelSize;
            this->stddev = stddev;
        }

        std::vector<Eigen::MatrixXd> apply();

};

class MedianFilter: public IKernel{

    private:

        static Eigen::MatrixXd filter(Eigen::MatrixXd A, int windowSize);

    public:

        MedianFilter(std::vector<Eigen::MatrixXd> channels) {this->channels = channels;}

        std::vector<Eigen::MatrixXd> applyFilter(std::vector<Eigen::MatrixXd> channels, int windowSize);

        std::vector<Eigen::MatrixXd> getChannels() {return this->channels;}

};