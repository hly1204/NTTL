#ifndef FP_HPP
#define FP_HPP

#include "field.hpp"
#include <cstdint>
#include <utility>
#include <type_traits>
#include <iostream>
#include <exception>
#include <limits>

namespace nttl {

/// \brief 域 \f$\mathbb{F}_p\f$ 其中 \f$p\f$ 为素数.
template<std::uint32_t P> requires is_prime_v<P> && (P >> 31 == 0)
class Fp {
public:
    struct UnsafeTag {};

    using value_type = std::uint32_t;
    static constexpr auto MOD = P;

    template<typename T>
    requires std::is_integral_v<T> && std::is_signed_v<T>
    static constexpr auto safe_mod(T v) {
        v %= static_cast<std::make_signed_t<std::decay_t<decltype(MOD)>>>(MOD);
        if (v < 0) v += static_cast<std::make_signed_t<std::decay_t<decltype(MOD)>>>(MOD);
        return v;
    }

    template<typename T>
    constexpr Fp(T v, UnsafeTag) requires std::is_integral_v<T>
        : v_(static_cast<std::decay_t<decltype(v_)>>(v)) {}
    constexpr Fp() {}
    template<typename T>
    constexpr Fp(T v) requires std::is_integral_v<T> && std::is_unsigned_v<T>
        : v_(static_cast<std::decay_t<decltype(v_)>>(v % MOD)) {}
    template<typename T>
    constexpr Fp(T v) requires std::is_integral_v<T> && std::is_signed_v<T>
        : v_(static_cast<std::decay_t<decltype(v_)>>(safe_mod(v))) {}
    static constexpr auto x() { return MOD; }
    static constexpr auto card() { return std::make_pair(MOD, 1u); }
    constexpr auto val() const { return v_; }

    /// \brief 加法逆元.
    constexpr auto operator-() const {
        if (v_ == 0) return Fp();
        return Fp(MOD - v_, UnsafeTag());
    }
    /// \brief 乘法逆元.
    constexpr auto inv() const {
        using S = std::make_signed_t<std::decay_t<decltype(v_)>>;
        if (v_ == 0) throw std::runtime_error("division by zero");
        S x1 = 1, x3 = 0, a = static_cast<S>(v_), b = static_cast<S>(MOD);
        while (b != 0) {
            S q = a / b, x1_old = x1, a_old = a;
            x1 = x3, x3 = x1_old - x3 * q, a = b, b = a_old - b * q;
        }
        return Fp(x1);
    }

    constexpr decltype(auto) operator+=(const Fp &rhs) {
        if ((v_ += rhs.v_) >= MOD) v_ -= MOD;
        return (*this);
    }
    constexpr decltype(auto) operator-=(const Fp &rhs) {
        if ((v_ = v_ + MOD - rhs.v_) >= MOD) v_ -= MOD;
        return (*this);
    }
    constexpr decltype(auto) operator*=(const Fp &rhs) {
        v_ = static_cast<std::uint64_t>(v_) * rhs.v_ % MOD;
        return (*this);
    }
    constexpr decltype(auto) operator/=(const Fp &rhs)
    try { return (*this *= rhs.inv()); } catch(...) { throw; }
    constexpr auto pow(long long e) const {
        if ((e %= MOD - 1) < 0) e += MOD - 1;
        for (Fp res{1}, x{*this};; x *= x) {
            if (e & 1) res *= x;
            if ((e >>= 1) == 0) return res;
        }
    }
    constexpr auto operator^(long long e) const { return pow(e); }
    constexpr decltype(auto) operator[](int n) {
        if (n != 0) throw std::out_of_range("only Fp<P>[0] will work");
        return (*this);
    }
    constexpr decltype(auto) operator[](int n) const {
        if (n != 0) throw std::out_of_range("only Fp<P>[0] will work");
        return (*this);
    }

    friend constexpr auto operator+(const Fp &lhs, const Fp &rhs) { return Fp(lhs) += rhs; }
    friend constexpr auto operator-(const Fp &lhs, const Fp &rhs) { return Fp(lhs) -= rhs; }
    friend constexpr auto operator*(const Fp &lhs, const Fp &rhs) { return Fp(lhs) *= rhs; }
    friend constexpr auto operator/(const Fp &lhs, const Fp &rhs)
    try { return Fp(lhs) /= rhs; } catch(...) { throw; }

    friend constexpr auto operator==(const Fp &lhs, const Fp &rhs) { return lhs.v_ == rhs.v_; }
    friend constexpr auto operator!=(const Fp &lhs, const Fp &rhs) { return lhs.v_ != rhs.v_; }

    friend decltype(auto) operator>>(std::istream &lhs, Fp &rhs) {
        using S = std::make_signed_t<std::decay_t<decltype(rhs.val())>>;
        S v;
        lhs >> v;
        rhs = Fp(v);
        return (lhs);
    }

    friend decltype(auto) operator<<(std::ostream &lhs, const Fp &rhs) {
        return (lhs << rhs.val());
    }

protected:
    std::uint32_t v_{};
};

}

#endif
