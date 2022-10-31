cmake_minimum_required(VERSION 3.12)
project(BalanceQP)
set(CMAKE_BUILD_TYPE "Debug")
# g++编译选项 
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -O3 -std=c++14 -Wall")

find_package(OsqpEigen) # 找OsqpEigen库
find_package(Eigen3) # 找Eigen3库

# 添加头文件路径 -I
include_directories(${CMAKE_SOURCE_DIR}/include
                    ${EIGEN3_INCLUDE_DIR}
)

# 添加源文件路径
set(SRC
    src/balance_qp_solver.cpp
) 

# 输出可执行文件 
add_executable(main src/check.cpp ${SRC})

# 链接库 -l
target_link_libraries(main PRIVATE 
                        pthread           # C线程库
                        OsqpEigen::OsqpEigen # OsqpEigen库
)
