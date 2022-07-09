cmake_minimum_required(VERSION 3.16)

# Number Theory Template Library
project(NTTL VERSION 0.0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (NOT MSVC)
    add_compile_options(-O0)
    add_compile_options(-Wall)
    add_compile_options(-Werror)
    add_compile_options(-fsanitize=undefined,address)
    add_link_options(-fsanitize=undefined,address)
else ()
    add_compile_options(/utf-8)
    add_compile_options(/Od)
    add_compile_options(/W4)
endif ()

add_subdirectory(test/unit_test)