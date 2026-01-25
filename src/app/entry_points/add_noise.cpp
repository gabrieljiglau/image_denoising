#include "../../libs/CLI11.hpp"
#include "../../libs/lodepng.h"
#include "../include/utils.hpp"
#include <string>

int main(int argc, char *argv[]){

    CLI::App app{"Add artificial noise (gaussian or salt-and-pepper) from the command line"};

    std::string inputPng;
    std::vector<std::string> outputPngs;
    
    app.add_option("--i, --input", inputPng, "Input image path (png format)")
        ->required();
    app.add_option("--o, --output", outputPngs, "Output paths (png format)")
        ->required();

    
    double mean = 0;
    std::vector<double> stddevs;
    auto *gaussianNoise = app.add_subcommand("gaussian_noise", "Gaussian Noise");
    gaussianNoise->add_option("--mean", mean, "The mean");
    gaussianNoise->add_option("--stddev", stddevs, "The standard deviation")
        ->required()
        ->check(CLI::PositiveNumber);

    double threshold = 0;
    auto *saltAndPepper = app.add_subcommand("salt_pepper", "Salt and Pepper Noise");
    saltAndPepper->add_option("--threshold", threshold, "The proportion of defunct pixels in [0, 1] interval")
        ->required()
        ->check(CLI::PositiveNumber);

    CLI11_PARSE(app, argc, argv);

    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;

    unsigned error = lodepng::decode(image, width, height, inputPng);

    if (*gaussianNoise){
        if (outputPngs.size() != stddevs.size()){
            std::cout << "There must be 1:1 ration between stddevs and the paths to store the image ! \n";
            return 1;
        }

        addGaussianNoise(image, height, width, mean, stddevs, outputPngs);

    } else if(*saltAndPepper){
        if (outputPngs.size() != 1){
            std::cout << "Only one image can be generated at once with salt and pepper noise ! \n";
            return 2;
        }

        if (threshold > 0.5) {
            std::cout << "We reccommend that threshold is inside [0, 0.5], otherwise the image is too messy ! \n";
            return 3;
        }
        
        addSaltPepperNoise(image, outputPngs[0], height, width, threshold);
    }

    return 0;
}