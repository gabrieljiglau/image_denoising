#include "../libs/CLI11.hpp"
#include "../include/toolkit.hpp"
#include <Eigen/src/Core/Matrix.h>
#include <Eigen/src/Core/util/Constants.h>
#include <vector>
#include <string>
#include <iostream>
#include <fmt/core.h>
#include <fmt/core.h>


int main(int argc, char *argv[]){

    CLI::App app{"Image denoising tool"};

    std::string inputPng;
    std::string outputPng;
    std::string mode;

    
    app.add_option("-i, --input", inputPng, "Input image path (png format)")
        ->required();
    app.add_option("--o, --output", outputPng, "Output path (png format)")
        ->required();
    app.add_option("--m, --mode", mode, "What algorithm to run")
        ->required();
    
    // Gaussian Blur mode
    int kernelSize;
    double stddev;

    auto *blur = app.add_subcommand("blur", "Gaussian Blur");
    blur->add_option("--kernelSize", kernelSize, "The size of the kernel")
        ->required()
        ->check(CLI::PositiveNumber);
    blur->add_option("--stddev", stddev, "The standard deviation")
        ->required();

    
    // Median Filter
    int windowSize;

    auto *filter = app.add_subcommand("filter", "Median Filter");
    filter->add_option("--windowSize", kernelSize, "The size of the kernel")
        ->required()
        ->check(CLI::PositiveNumber);

    
    // SVD mode
    std::vector<int> kVals;
    int maxIterations; // 10000
    double epsilon; // 1e-10

    auto *svd = app.add_subcommand("svd", "Singular Value Decomposition");
    svd->add_option("--k", kVals, "Keep only the k-biggest singular values")
        ->required()
        ->check(CLI::PositiveNumber);
    svd->add_option("--maxIterations", maxIterations, "Maximum Iterations")
        ->required()
        ->check(CLI::PositiveNumber);
    svd->add_option("--epsilon", epsilon, "Precision after the decimal point")
        ->required();

    // SVD with patches
    int patchSize;
    int stride;

    /// TODO: sa vezi ca merge cmake-ul pentru main-ul algoritmului
    /// TODO: apoi sa dai ca vector (mai multe valori pentru optiuni) pt argumentele blur si filter

    auto *patchSvd = app.add_subcommand("patch_svd", "Singular Value Decomposition with patches");
    patchSvd->add_option("--k", kVals, "Keep only the k-biggest singular values")
        ->required()
        ->check(CLI::PositiveNumber);
    patchSvd->add_option("--maxIterations", maxIterations, "Maximum Iterations")
        ->required()
        ->check(CLI::PositiveNumber);
    patchSvd->add_option("--epsilon", epsilon, "Precision after the decimal point")
        ->required();
    patchSvd->add_option("--patchSize", patchSize, "The dimensions of the smaller submatrices (patchSize X patchSize)t")
        ->required();
    patchSvd->add_option("--stride", stride, "The step size then traversing the channel")
        ->required();
        
    CLI11_PARSE(app, argc, argv);

    Toolkit toolkit(inputPng, mode);
    
    if (*blur) {
        toolkit.setKernelSize(kernelSize);
        toolkit.setStddev(stddev);
        
        std::cout << "Applying Gaussian Blur \n" << std::endl;
    
    } else if (*filter){
        toolkit.setWindowSize(windowSize);
        
        std::cout << "Applying Median Filter \n" << std::endl;
    
    } else if(*svd || *patchSvd){
        toolkit.setMaxIterations(maxIterations);
        toolkit.setEpsilon(epsilon);
        toolkit.checkRank();

        if (*patchSvd){
            toolkit.setPatchSize(patchSize);
            toolkit.setStride(stride);

            std::cout << "Applying patch SVD \n" << std::endl;
        } else {
            std::cout << "Applying standard SVD \n" << std::endl;
        }
    }

    std::vector<int> errCodes;
    if (!*svd || *patchSvd){
        errCodes = toolkit.processPng(inputPng, outputPng);   
    } else {
        toolkit.setK(kVals);
        errCodes = toolkit.processPng(inputPng, outputPng);
    }

    for (int errCode : errCodes){
        if (errCode != 0) {
            std::cout << "Encountered error code : " << errCode << " when processing the image with lodepng" << std::endl;
        }
    }

    return 0;

    /// TODO: apoi rezolvat cmake-ul
    /// TODO: de scris niste teste cu catch 2

    /*
    std::vector<int> stdDevs = {10};
    addGaussianNoise(image, height, width, stdDevs);
    */
    /*
    addSaltPepperNoise(image, newImage, "/home/gabriel/Documents/HolyC/image_denoising/images/noisy/car_salt_pepper.png",
                     height, width, 0.1);
    */
    /*
    std::vector<Eigen::MatrixXd> channels = rgbChannel(image, height, width);
    int windowSize = 3;
    std::vector<Eigen::MatrixXd> medianChannels = applyMedianFilter(channels, windowSize);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_car_median.png");
    reconstructImage(medianChannels, newImage, newPath, height, width);
    */  
    
    
    // gaussianBlur after SVD
    /*
    int kernelSize = 3;
    std::vector<int> stdDevs = {2, 5, 10, 20};
    for (int stdDev : stdDevs){
        std::vector<Eigen::MatrixXd> blurredChannels = applyGaussianBlur(channels, kernelSize, stdDev);
        std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised/denoised_k{}_afterBlur.png", stdDev);
        reconstructImage(blurredChannels, newImage, newPath, height, width);
    }
    */



    /*
    // median filter (before/after SVD); still bad results
    int windowSize = 2;
    std::vector<Eigen::MatrixXd> medianChannels = applyMedianFilter(channels, windowSize);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k50_median.png");
    reconstructImage(medianChannels, newImage, newPath, height, width);
    */

    // svd on smaller submatrices of the original image, bad results
    /*=
    int patchSize = 16;
    int stride = 12;
    int k = 2;
    
    std::vector<Eigen::MatrixXd> finalChannels = patchesDecomp(channels, patchSize, stride, k, maxIterations, epsilon);
    std::string newPath = fmt::format("/home/gabriel/Documents/HolyC/image_denoising/images/denoised_k{}.png", k);
    reconstructImage(finalChannels, newImage, newPath, height, width);
    */

}