#ifndef FQ_HPP
#define FQ_HPP

#include "field.hpp"
#include "polynomial.hpp"
#include "fp.hpp"
#include <array>
#include <cstdint>
#include <utility>
#include <exception>

namespace nttl {

/// \brief 域 \f$\mathbb{F}_q\f$ 其中 \f$q\f$ 为素数的幂次.
template<FiniteField FF, std::size_t DEG, const std::array<FF, DEG + 1> &IrrPoly>
    requires (FF::card().second == 1)
class Fq : public std::array<FF, DEG> {
public:
    using MyBase = std::array<FF, DEG>;
    using MyBase::array;
    using value_type = typename MyBase::value_type;

    constexpr Fq(const value_type &v) { MyBase::front() = v; }
    constexpr Fq(const Fq &) = default;
    constexpr Fq &operator=(const Fq &) = default;

    operator Poly<FF>() const { return Poly(MyBase::cbegin(), MyBase::cend()); }

    static constexpr auto SIZE = static_cast<int>(DEG);
    static constexpr auto POLY = IrrPoly;

    static constexpr auto x() { return value_type::x(); }
    static constexpr auto card() { return std::make_pair(x(), static_cast<std::uint32_t>(DEG)); }

    enum : int {
        NEGATIVE_INFINITY = -1,
    };

    constexpr auto deg() const {
        auto d = SIZE - 1;
        while (d >= 0 && MyBase::operator[](d) == value_type{}) --d;
        return d < 0 ? static_cast<int>(NEGATIVE_INFINITY) : d;
    }
    /// \brief 首项系数
    constexpr auto lc() const {
        const auto d = deg();
        return d == NEGATIVE_INFINITY ? value_type{} : MyBase::operator[](d);
    }

    /// \brief 加法逆元.
    constexpr auto operator-() const {
        Fq res(*this);
        for (auto &&i : res) i = -i;
        return res;
    }
    /// \brief 乘法逆元
    constexpr auto inv() const {
        using P = Poly<value_type>;
        if (*this == Fq{}) throw std::runtime_error("Division by zero");
        P iv;
        std::tie(iv, std::ignore) = P::inv_gcd(*this, P(POLY.cbegin(), POLY.cend()));
        Fq res;
        for (auto i = 0, e = iv.deg(); i <= e; ++i) res[i] = iv[i];
        return res;
    }
    constexpr decltype(auto) operator+=(const Fq &rhs) {
        for (auto i = 0; i != SIZE; ++i) MyBase::operator[](i) += rhs[i];
        return (*this);
    }
    constexpr decltype(auto) operator-=(const Fq &rhs) {
        for (auto i = 0; i != SIZE; ++i) MyBase::operator[](i) -= rhs[i];
        return (*this);
    }
    constexpr decltype(auto) operator*=(const Fq &rhs) {
        using P = Poly<value_type>;
        P res(P(*this) * P(rhs) % P(POLY.cbegin(), POLY.cend()));
        const auto d = res.deg();
        for (auto i = 0; i <= d; ++i) MyBase::operator[](i) = res[i];
        for (auto i = d + 1; i < SIZE; ++i) MyBase::operator[](i) = value_type{}; // 重要!!!
        return (*this);
    }
    constexpr decltype(auto) operator/=(const Fq &rhs)
    try { return (*this *= rhs.inv()); } catch(...) { throw; }
    constexpr auto pow(long long e) const {
        Fq res{1}, x(*this);
        if (e < 0) {
            x = x.inv();
            e = -e; // 为了简便, 不考虑溢出
        }
        for (;;) {
            if (e & 1) res *= x;
            if ((e >>= 1) == 0) return res;
            x *= x;
        }
    }
    constexpr auto operator^(long long e) const { return pow(e); }

    friend constexpr auto operator+(const Fq &lhs, const Fq &rhs) { return Fq(lhs) += rhs; }
    friend constexpr auto operator-(const Fq &lhs, const Fq &rhs) { return Fq(lhs) -= rhs; }
    friend constexpr auto operator*(const Fq &lhs, const Fq &rhs) { return Fq(lhs) *= rhs; }
    friend constexpr auto operator/(const Fq &lhs, const Fq &rhs)
    try { return Fq(lhs) /= rhs; } catch(...) { throw; }

    friend decltype(auto) operator<<(std::ostream &lhs, const Fq &rhs) {
        auto s = 0;
        lhs << '[';
        for (auto &&i : rhs) {
            lhs << i;
            if (s >= 1) lhs << 'x';
            if (s > 1) lhs << '^' << s;
            if (++s != SIZE) lhs << " + ";
        }
        return (lhs << ']');
    }
};

namespace detail {

/// \see https://math.stackexchange.com/questions/4092518/construct-irreducible-polynomials-of-degree-32-over-z-2x
static constexpr std::array<Fp<2>, 9> IrrPoly_2_8{1, 1, 0, 1, 1, 0, 0, 0, 1};
static constexpr std::array<Fp<2>, 33> IrrPoly_2_32{1, 0, 0, 1, 1, 0, 0, 1,
                                                    0, 1, 0, 0, 0, 0, 0, 1,
                                                    0, 0, 0, 0, 0, 0, 0, 0,
                                                    0, 0, 0, 0, 0, 0, 0, 0, 1};

}

using F_2_8 = Fq<Fp<2>, 8, detail::IrrPoly_2_8>;
using F_2_32 = Fq<Fp<2>, 32, detail::IrrPoly_2_32>;

}

#endif
