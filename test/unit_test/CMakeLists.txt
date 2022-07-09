cmake_policy(SET CMP0135 NEW)
# see https://google.github.io/googletest/quickstart-cmake.html
include(FetchContent)
FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)
# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

# 该宏定义需要到对应路径下的调用 `ctest` 命令, 而不是在构建的目录下.
enable_testing()

include(GoogleTest)

file(GLOB test_sources *.cpp)

foreach (X IN LISTS test_sources)
    message("Add to testcases: " ${X})
    get_filename_component(FILENAME_WITHOUT_EXT ${X} NAME_WE)
    add_executable(${FILENAME_WITHOUT_EXT} ${X})
    target_include_directories(${FILENAME_WITHOUT_EXT} PUBLIC ${CMAKE_SOURCE_DIR}/src)
    target_link_libraries(${FILENAME_WITHOUT_EXT} PRIVATE gtest_main)
    gtest_discover_tests(${FILENAME_WITHOUT_EXT} DISCOVERY_TIMEOUT 120)
endforeach ()
