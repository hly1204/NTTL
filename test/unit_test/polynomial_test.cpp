#include "fq.hpp"
#include "polynomial.hpp"
#include "random.hpp"
#include <gtest/gtest.h>
#include <random>

TEST(Poly, BasicTest) {
    nttl::xoshiro256starstar gen(std::random_device{}());
    auto r = [&gen]() {
        nttl::F_2_32 a;
        std::uniform_int_distribution<int> dis(0, 1);
        for (auto i = 0; i != static_cast<int>(a.size()); ++i) a[i] = dis(gen);
        return a;
    };
    nttl::Poly<nttl::F_2_32> f, g;
    for (auto i = 0; i != 10; ++i) f.emplace_back(r());
    for (auto i = 0; i != 10; ++i) g.emplace_back(r());
    auto k = r();
    EXPECT_EQ(f(k) * g(k), (f * g)(k)) << "f(k) * g(k) == (fg)(k)";
}

TEST(Poly, InterpolationTest1) {
    nttl::xoshiro256starstar gen(std::random_device{}());
    auto r = [&gen]() {
        nttl::F_2_32 a;
        std::uniform_int_distribution<int> dis(0, 1);
        for (auto i = 0; i != static_cast<int>(a.size()); ++i) a[i] = dis(gen);
        return a;
    };
    nttl::Poly<nttl::F_2_32> f;
    // deg(`f`) = 9
    for (auto i = 0; i != 10; ++i) f.emplace_back(r());
    std::vector<nttl::F_2_32> x, y;
    for (auto i = 0; i != 10; ++i) {
        nttl::F_2_32 v;
        for (auto j = 0; j != static_cast<int>(v.size()); ++j) v[j] = i >> j & 1;
        x.emplace_back(v);
        y.emplace_back(f(v));
    }
    auto [ff, m] = f.inter(x, y);
    EXPECT_EQ(ff, f) << "f == ff";
}

TEST(Poly, InterpolationTest2) {
    nttl::xoshiro256starstar gen(std::random_device{}());
    std::uniform_int_distribution<int> dis(0, 998244352);
    static constexpr std::array<nttl::Fp<998244353>, 2> IrrP{0, 1};
    nttl::Poly<nttl::Fq<nttl::Fp<998244353>, 1, IrrP>> f;
    // deg(`f`) = 9
    for (auto i = 0; i != 10; ++i) f.emplace_back(dis(gen));
    std::vector<nttl::Fq<nttl::Fp<998244353>, 1, IrrP>> x, y;
    for (auto i = 0; i != 10; ++i) {
        nttl::Fq<nttl::Fp<998244353>, 1, IrrP> v{i};
        x.emplace_back(v);
        y.emplace_back(f(v));
    }
    auto [ff, m] = f.inter(x, y);
    EXPECT_EQ(ff, f) << "f == ff";
}

TEST(Poly, InterpolationWithErrorTest) {
    nttl::xoshiro256starstar gen(std::random_device{}());
    auto r = [&gen]() {
        nttl::F_2_32 a;
        std::uniform_int_distribution<int> dis(0, 1);
        for (auto i = 0; i != static_cast<int>(a.size()); ++i) a[i] = dis(gen);
        return a;
    };
    nttl::Poly<nttl::F_2_32> f;
    const auto kp = 10; // deg(`f`) = 9
    const auto l = 3, k = 2 * l + kp;
    for (auto i = 0; i != kp; ++i) f.emplace_back(r());
    std::vector<nttl::F_2_32> x, y;
    for (auto i = 0; i != k; ++i) {
        nttl::F_2_32 v;
        for (auto j = 0; j != static_cast<int>(v.size()); ++j) v[j] = i >> j & 1;
        x.emplace_back(v);
        y.emplace_back(f(v));
    }
    for (int i = 0; i != 10; ++i) {
        std::uniform_int_distribution<int> dis(0, k - 1);
        auto yy = y;
        for (auto j = 0; j != l; ++j) yy[dis(gen)] = r();
        EXPECT_EQ(f, f.inter_we(x, yy, kp, l).value_or(nttl::Poly<nttl::F_2_32>{}));
    }
}

int main(int argc, char *argv[]) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}