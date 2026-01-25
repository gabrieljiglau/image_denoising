#include "catch2/catch_test_macros.hpp"
#include "../src/app/include/toolkit.hpp"
#include "../src/libs/lodepng.h"
#include "../src/app/include/utils.hpp"


TEST_CASE("Input that exists, with a different extension than PNG", "[utils]"){

    bool checkExistence = true;
    int returnCode = validateInput("images/lena.jpg", checkExistence);
    REQUIRE(returnCode == 1);
}

TEST_CASE("Input of the same type doesn't exist", "[utils]"){

    bool checkExistence = true;
    int returnCode = validateInput("images/404_error.png", checkExistence);
    REQUIRE(returnCode == 2);
}

TEST_CASE("Input of different type doesn't exist", "[utils]"){

    bool checkExistence = true;
    int returnCode = validateInput("images/404_error.jpeg", checkExistence);
    REQUIRE(returnCode == 1);
}


TEST_CASE("Big Integer for stddev", "[gaussian_noise]"){

    std::string inputPng = "images/lena.png";
    std::vector<std::string> outputPngs = {"images/lena_big_stddev.png"};
    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;
    
    double mean = 0;
    std::vector<double> stddevs = {500000};

    unsigned error = lodepng::decode(image, width, height, inputPng);
    int returnCode = addGaussianNoise(image, height, width, mean, stddevs, outputPngs);

    // image looks demonic, but that's life

    REQUIRE(returnCode == 0);
}

TEST_CASE("Big Integer for mean", "[gaussian_noise]"){

    std::string inputPng = "images/lena.png";
    std::vector<std::string> outputPngs = {"images/lena_big_mean.png"};
    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;
    
    double mean = 50000;
    std::vector<double> stddevs = {0};

    unsigned error = lodepng::decode(image, width, height, inputPng);
    int returnCode = addGaussianNoise(image, height, width, mean, stddevs, outputPngs);

    REQUIRE(returnCode == 0);
}

TEST_CASE("Multiple stddevs", "[gaussian_noise]"){

    std::string inputPng = "images/lena.png";
    std::vector<std::string> outputPngs = {"images/lena_1.png", "images/lena_2.png", "images/lena_3.png"};
    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;
    
    double mean = -2;
    std::vector<double> stddevs = {-3, 3, 6};

    unsigned error = lodepng::decode(image, width, height, inputPng);
    int returnCode = addGaussianNoise(image, height, width, mean, stddevs, outputPngs);

    REQUIRE(returnCode == 0);
}


TEST_CASE("Smallest threshold", "[salt_pepper_noise]"){

    std::string inputPng = "images/lena.png";
    std::vector<std::string> outputPngs = {"images/lena_small_threshold.png"};
    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;
    
    double threshold = 0;
    // it works, but it's undesirable; the final image will be meaningless, made up of only white / black pixels

    unsigned error = lodepng::decode(image, width, height, inputPng);
    int returnCode = addSaltPepperNoise(image, outputPngs[0], height, width, threshold);

    REQUIRE(returnCode == 0);
}

TEST_CASE("Out of bounds threshold", "[salt_pepper_noise]"){

    std::string inputPng = "images/lena.png";
    std::vector<std::string> outputPngs = {"images/lena_big_threshold.png"};
    std::vector<unsigned char> image;
    unsigned int height;
    unsigned int width;
    
    double threshold = 0;

    unsigned error = lodepng::decode(image, width, height, inputPng);
    int returnCode = addSaltPepperNoise(image, outputPngs[0], height, width, threshold);

    REQUIRE(returnCode == 0);
}