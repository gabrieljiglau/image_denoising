#pragma once
#include <optional>
#include <string>
#include <vector>
#include "lodepng.h"

class Toolkit{

    private:
        std::vector<unsigned char> image;
        int width;
        int height;
        std::string mode;
        std::optional<std::string> kSVD;

        int loadPng(std::string pngPath){
            // out image for reconstruction
            std::vector<unsigned char> newImage;
            unsigned width; // img_row
            unsigned height; // img_col

            unsigned error = lodepng::decode(image, width, height, pngPath);
            return error;
        }

    public:
        Toolkit(){
            this->width = 0;
            this->height = 0;
        }

        

}