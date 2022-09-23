cmake_minimum_required(VERSION 3.12)
project(QP_SOLVE)
set(CMAKE_BUILD_TYPE "Debug")
# g++编译选项 
set(CMAKE_CXX_FLAGS "$ENV{CXXFLAGS} -O3 -std=c++14 -Wall")

# using pkg-config find pinocchio
set(ENV{PKG_CONFIG_PATH} /opt/openrobots/lib/pkgconfig)
find_package(PkgConfig REQUIRED) 
pkg_check_modules(PIN REQUIRED pinocchio) # PIN 
message(STATUS "=== PIN_LIBRARIES: ${PIN_LIBRARIES}") # lib info
message(STATUS "=== PIN_INCLUDE_DIRS: ${PIN_INCLUDE_DIRS}") # pinocchio include path

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

# -L 
link_directories(/opt/openrobots/lib) # .so file path

# 输出可执行文件 
add_executable(main src/QP_solve.cpp ${SRC})

# 链接库 -l
target_link_libraries(main PRIVATE 
                        pthread           # C线程库
                        OsqpEigen::OsqpEigen # OsqpEigen库
)

target_include_directories(main PRIVATE ${PIN_INCLUDE_DIRS})
target_compile_options(main PRIVATE ${PIN_CFLAGS_OTHER})
