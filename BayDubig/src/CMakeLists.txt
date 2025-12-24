cmake_minimum_required(VERSION 3.26)
project(src)

set(CMAKE_CXX_STANDARD 11)  # Change to 14 if needed
add_executable(src Baydubig_mcmc.cpp)

set(ARMADILLO_INCLUDE_DIR /opt/homebrew/Cellar/armadillo/12.6.6/include)
set(ARMADILLO_LIBRARY /opt/homebrew/lib/libarmadillo.dylib)

# Set the R installation prefix
set(CMAKE_PREFIX_PATH "/usr/lib/R")  # Ensure this is the correct path
set(R_DIR "/usr/lib/R")  # Ensure this is the correct path

# Find packages
find_package(Armadillo REQUIRED)
find_package(R REQUIRED)
find_package(Rcpp REQUIRED)
find_package(RcppArmadillo REQUIRED)
find_package(OpenMP REQUIRED)
target_link_libraries(src PUBLIC OpenMP::OpenMP_CXX)

# Include directories
include_directories(${ARMADILLO_INCLUDE_DIR})
include_directories(${Rcpp_INCLUDE_DIRS})
include_directories(${RcppArmadillo_INCLUDE_DIRS})
include_directories("/usr/share/R/include")  # R header files

# Link libraries
target_link_libraries(src ${Armadillo_LIBRARY} ${Rcpp_LIBRARIES} ${RcppArmadillo_LIBRARIES})
