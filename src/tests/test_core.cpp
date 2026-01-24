#include "catch2/catch_test_macros.hpp"
#include "../app/include/toolkit.hpp"
#include "../libs/lodepng.h"

/// fails, however, added safeguards in main; the program returns if the .png extension is not present
TEST_CASE("Input with a different extension than PNG", "[tookit]"){

    Toolkit toolkit = Toolkit("images/lena.jpeg");
    toolkit.initChannels();
}