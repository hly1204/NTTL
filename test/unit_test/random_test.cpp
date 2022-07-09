#include "fq.hpp"
#include "polynomial.hpp"
#include "random.hpp"
#include <gtest/gtest.h>
#include <random>

TEST(xoshiro256starstar, BasicTest) {
    nttl::xoshiro256starstar gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, 10);
    for (int i = 0; i != 10; ++i) {
        int r = dis(gen);
        EXPECT_LE(r, 10);
        EXPECT_GE(r, 0);
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}