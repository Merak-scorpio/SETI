cmake_minimum_required(VERSION 3.10)

project(SETI)

include_directories(include)
aux_source_directory(./src SRC_FILES)

# 设置指定的C++编译器版本是必须的，如果不设置，或者为OFF，则指定版本不可用时，会使用上一版本。
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 指定为C++标准 版本
set(CMAKE_CXX_STANDARD 20)

#设定编译后生成的可执行二进制文件所在的目录
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)

#设定编译后生成的库文件所在的目录
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
#并且把该目录设为链接目录
link_directories(${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp -Wall")

set(CMAKE_BUILD_TYPE Release)

add_executable(SETI main.cpp ${SRC_FILES})