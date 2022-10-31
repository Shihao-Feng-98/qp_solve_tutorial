cmake_minimum_required(VERSION 3.12)
project(ConvexMPC)

# g++编译选项 -O3优化代码
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -O3 -std=c++14 -Wall")

# spdlog 
if(NOT TARGET spdlog)
    # Stand-alone build
    find_package(spdlog REQUIRED)
endif()
find_package(OsqpEigen REQUIRED) 
find_package(Eigen3) 

# 添加头文件路径 -I
include_directories(${CMAKE_SOURCE_DIR}/include    
                    ${EIGEN3_INCLUDE_DIR}
)

# 添加源文件路径
set(SRC 
    src/ConvexMpcSolver.cpp
)

# 输出可执行文件
add_executable(main src/check.cpp ${SRC})

# 链接库 -l
target_link_libraries(main PRIVATE 
                    pthread           # C线程库
                    OsqpEigen::OsqpEigen # OsqpEigen lib
                    spdlog::spdlog       # spdlog lib
)
