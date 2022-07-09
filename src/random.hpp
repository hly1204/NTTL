#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <cstdint>
#include <limits>

namespace nttl {

/// \brief 随机数类.
/// \see https://prng.di.unimi.it/xoshiro256starstar.c
/// \remarks original license CC0 1.0
class xoshiro256starstar {
public:
    /// \see https://en.cppreference.com/w/cpp/named_req/UniformRandomBitGenerator
    using result_type = std::uint64_t;

    /// \see https://prng.di.unimi.it/splitmix64.c
    /// \remarks original license CC0 1.0
    explicit xoshiro256starstar(result_type seed) {
        for (int i = 0; i != 4; ++i) {
            result_type z = (seed += 0x9e3779b97f4a7c15);
            z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9;
            z = (z ^ (z >> 27)) * 0x94d049bb133111eb;
            s_[i] = z ^ (z >> 31);
        }
    }
    static constexpr auto min() { return std::numeric_limits<result_type>::min(); }
    static constexpr auto max() { return std::numeric_limits<result_type>::max(); }

private:
    static inline auto rotl(const result_type x, int k) { return (x << k) | (x >> (64 - k)); }
    auto next() {
        const auto res = rotl(s_[1] * 5, 7) * 9;
        const auto t = s_[1] << 17;
        s_[2] ^= s_[0];
        s_[3] ^= s_[1];
        s_[1] ^= s_[2];
        s_[0] ^= s_[3];
        s_[2] ^= t;
        s_[3] = rotl(s_[3], 45);
        return res;
    }
    result_type s_[4];

public:
    auto operator()() { return next(); }
};

}

#endif
