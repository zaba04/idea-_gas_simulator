cmake_minimum_required(VERSION 3.30)
project(kms_1)

set(CMAKE_CXX_STANDARD 20)

add_executable(kms_1 main.cpp
        atom.cpp
        atom.h
        parameters.cpp
        parameters.h)
