#include "../../libs/CLI11.hpp"
#include "../include/toolkit.hpp"
#include "../include/utils.hpp"
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
    std::vector<std::string> outputPngs;
    
    app.add_option("-i, --input", inputPng, "Input image path (png format)")
        ->required();
    app.add_option("--o, --output", outputPngs, "Output paths (png format)")
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
    auto *filter = app.add_subcommand("median_filter", "Median Filter");
    filter->add_option("--kernelSize", kernelSize, "The size of the kernel")
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

    // validate input-output
    bool checkExistence = true;
    validateInput(inputPng, checkExistence); /*inputPng, checkExistence*/

    for (std::string &path : outputPngs){
        validateInput(path, !checkExistence);
    }


    Toolkit toolkit(inputPng);
    toolkit.initChannels();
    
    if (*blur) {
        toolkit.setKernelSize(kernelSize);
        toolkit.setStddev(stddev);
        toolkit.setMode("blur");
        
        std::cout << "Applying Gaussian Blur \n" << std::endl;
    
    } else if (*filter){
        toolkit.setWindowSize(kernelSize);
        toolkit.setMode("median_filter");
        
        std::cout << "Applying Median Filter \n" << std::endl;
    
    } else if(*svd || *patchSvd) {
        toolkit.setK(kVals);
        toolkit.setMaxIterations(maxIterations);
        toolkit.setEpsilon(epsilon);
        toolkit.checkRank();

        if (*patchSvd){
            toolkit.setPatchSize(patchSize);
            toolkit.setStride(stride);
            toolkit.setMode("patch_svd");

            std::cout << "Applying patch SVD \n" << std::endl;
        } else {
            std::cout << "Applying standard SVD \n" << std::endl;
            toolkit.setMode("svd");
        }
    }

    std::vector<int> errCodes;
    errCodes = toolkit.processPng(inputPng, outputPngs);   

    for (int errCode : errCodes){
        if (errCode != 0) {
            std::cout << "Encountered error code : " << errCode << " when processing the image with lodepng" << std::endl;
        }
    } 

    return 0;

}