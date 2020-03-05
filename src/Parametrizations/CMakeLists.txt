set(CMAKE_CXX_FLAGS "-std=c++11 -O3 -Wall -Wextra")

include_directories(${PROJECT_SOURCE_DIR}/include)

FILE(GLOB ParamSources ./*.cpp)

add_library(parametrizations STATIC ${ParamSources})

