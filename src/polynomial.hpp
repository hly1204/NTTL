#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include "field.hpp"
#include <vector>
#include <exception>
#include <utility>
#include <tuple>
#include <iterator>
#include <iostream>
#include <optional>

namespace nttl {

/// \brief 多项式.
template<Field F>
class Poly : public std::vector<F> {
public:
    using MyBase = std::vector<F>;
    using MyBase::vector;
    using value_type = typename MyBase::value_type;

    enum : int {
        NEGATIVE_INFINITY = -1,
    };

    constexpr Poly(const Poly &) = default;
    constexpr Poly &operator=(const Poly &) = default;

    constexpr auto deg() const {
        auto d = static_cast<int>(MyBase::size()) - 1;
        while (d >= 0 && MyBase::operator[](d) == value_type{}) --d;
        return d < 0 ? static_cast<int>(NEGATIVE_INFINITY) : d;
    }
    constexpr decltype(auto) shrink() {
        MyBase::resize(deg() + 1);
        return (*this);
    }
    /// \brief 首项系数.
    constexpr auto lc() const {
        const auto d = deg();
        return d == NEGATIVE_INFINITY ? value_type{} : MyBase::operator[](d);
    }
    /// \brief 形式导数.
    /// \remarks 若在有限域中我们不用极限定义导数, 所以是形式的.
    constexpr auto deriv() const {
        const auto n = static_cast<int>(MyBase::size());
        if (n <= 1) return Poly{};
        Poly res;
        res.reserve(n - 1);
        for (auto i = 1; i != n; ++i) res.emplace_back(MyBase::operator[](i) * i);
        return res.shrink();
    }
    /// \brief 积分.
    constexpr auto integr(value_type c = value_type{}) const {
        const auto n = static_cast<int>(MyBase::size()) + 1;
        Poly res{c};
        res.reserve(n);
        for (auto i = 1; i != n; ++i) res.emplace_back(MyBase::operator[](i - 1) / i);
        return res.shrink();
    }
    constexpr auto operator-() const {
        Poly res;
        res.reserve(MyBase::size());
        for (auto &&i : *this) res.emplace_back(-i);
        return res.shrink();
    }
    constexpr decltype(auto) operator+=(const Poly &rhs) {
        if (MyBase::size() < rhs.size()) MyBase::resize(rhs.size());
        for (auto i = 0, e = static_cast<int>(rhs.size()); i != e; ++i) MyBase::operator[](i) += rhs[i];
        return (shrink());
    }
    constexpr decltype(auto) operator-=(const Poly &rhs) {
        if (MyBase::size() < rhs.size()) MyBase::resize(rhs.size());
        for (auto i = 0, e = static_cast<int>(rhs.size()); i != e; ++i) MyBase::operator[](i) -= rhs[i];
        return (shrink());
    }
    constexpr decltype(auto) operator*=(const Poly &rhs) {
        const auto n = deg(), m = rhs.deg();
        if (n == NEGATIVE_INFINITY || m == NEGATIVE_INFINITY) return (operator=(Poly{}));
        Poly res(n + m + 1);
        for (auto i = 0; i <= n; ++i)
            for (auto j = 0; j <= m; ++j)
                res[i + j] += MyBase::operator[](i) * rhs[j];
        return (operator=(res.shrink()));
    }
    constexpr decltype(auto) operator/=(const Poly &rhs) {
        auto n = deg(), m = rhs.deg(), q = n - m;
        if (m == -1) throw std::runtime_error("Division by zero");
        if (q <= -1) return (operator=(Poly{}));
        Poly res(q + 1);
        const auto iv = rhs.lc().inv();
        for (auto i = q; i >= 0; --i)
            if ((res[i] = MyBase::operator[](n--) * iv) != value_type{})
                for (auto j = 0; j != m; ++j) MyBase::operator[](i + j) -= res[i] * rhs[j];
        return (operator=(res));
    }
    constexpr decltype(auto) operator%=(const Poly &rhs) {
        auto n = deg(), m = rhs.deg(), q = n - m;
        if (m == -1) throw std::runtime_error("Division by zero");
        const auto iv = rhs.lc().inv();
        for (auto i = q; i >= 0; --i)
            if (auto res = MyBase::operator[](n--) * iv; res != value_type{})
                for (auto j = 0; j <= m; ++j) MyBase::operator[](i + j) -= res * rhs[j];
        return (shrink());
    }

    constexpr auto div_mod(const Poly &rhs) const {
        auto n = deg(), m = rhs.deg(), q = n - m;
        if (m == -1) throw std::runtime_error("Division by zero");
        if (q <= -1) return std::make_pair(Poly{}, Poly(*this).shrink());
        const auto iv = rhs.lc().inv();
        Poly quo(q + 1), rem(*this);
        for (auto i = q; i >= 0; --i)
            if ((quo[i] = rem[n--] * iv) != value_type{})
                for (auto j = 0; j <= m; ++j) rem[i + j] -= quo[i] * rhs[j];
        return std::make_pair(quo, rem.shrink()); // (quotient, remainder)
    }

    constexpr auto operator()(const value_type &pt) const {
        value_type res;
        for (auto i = deg(); i >= 0; --i) res = pt * res + MyBase::operator[](i);
        return res;
    }

    static constexpr auto inv_gcd(Poly a, Poly b) {
        Poly x1{value_type{1}}, x3;
        while (b != Poly{}) {
            auto [q, r] = a.div_mod(b);
            std::tie(x1, x3, a, b) = std::make_tuple(x3, x1 - x3 * q, b, r);
        }
        // 使 gcd 为首一多项式, 在 gcd 为一时可以直接获取乘法逆元
        return std::make_pair(x1 / Poly{a.lc()}, a / Poly{a.lc()});
    }

    /// \brief 插值.
    /// \param x 设多项式为 \f$f\f$ 那么 \c x[i] 为 \f$f\f$ 所求值的点.
    /// \param y 设多项式为 \f$f\f$ 那么 \c y[i] 为 \f$f\f$ 在 \c x[i] 的点值.
    static constexpr auto inter(const std::vector<value_type> &x, const std::vector<value_type> &y) {
        // `x` 中的元素必须唯一
        if (x.size() != y.size()) throw std::runtime_error("Size error");
        const auto n = static_cast<int>(x.size());
        Poly f, m{value_type{1}};
        for (auto i = 0; i != n; ++i) {
            f += Poly{(y[i] - f(x[i])) / m(x[i])} * m;
            m *= Poly{-x[i], value_type{1}};
        }
        return std::make_pair(f, m);
    }

    /// \brief 带有错误的插值.
    /// \param x  设多项式为 \f$f\f$ 那么 \c x[i] 为 \f$f\f$ 所求值的点. 共 \f$k\f$ 个.
    /// \param y  设多项式为 \f$f\f$ 那么 \c y[i] 为 \f$f\f$ 在 \c x[i] 的点值. 共 \f$k\f$ 个.
    /// \param kp \f$\deg(f) < k'\f$.
    /// \param l  最多有 \f$l\f$ 个位置的点值是错误的.
    static constexpr auto inter_we(const std::vector<value_type> &x, const std::vector<value_type> &y, int kp, int l)
        -> std::optional<Poly> {
        const auto k = static_cast<int>(x.size());
        if (k < 2 * l + kp) return {};
        auto [ff, m] = inter(x, y);
        auto extended_euclidean = [](Poly a, Poly b, int kp, int l) -> std::optional<Poly> {
            Poly x1{value_type{1}}, x2, x3, x4{value_type{1}};
            while (b != Poly{}) {
                auto [q, r] = a.div_mod(b);
                std::tie(x1, x2, x3, x4, a, b) = std::make_tuple(x3, x4, x1 - x3 * q, x2 - x4 * q, b, r);
                if (a.deg() - x2.deg() < kp && x2.deg() <= l && a % x2 == Poly{}) return a / x2;
            }
            return {};
        };
        return extended_euclidean(m, ff, kp, l);
    }

    friend constexpr auto operator+(const Poly &lhs, const Poly &rhs) { return Poly(lhs) += rhs; }
    friend constexpr auto operator-(const Poly &lhs, const Poly &rhs) { return Poly(lhs) -= rhs; }
    friend constexpr auto operator*(const Poly &lhs, const Poly &rhs) { return Poly(lhs) *= rhs; }
    friend constexpr auto operator/(const Poly &lhs, const Poly &rhs)
    try { return Poly(lhs) /= rhs; } catch(...) { throw; }
    friend constexpr auto operator%(const Poly &lhs, const Poly &rhs)
    try { return Poly(lhs) %= rhs; } catch(...) { throw; }

    friend constexpr auto operator==(const Poly &lhs, const Poly &rhs) {
        auto d = lhs.deg();
        if (d != rhs.deg()) return false;
        for (; d >= 0; --d)
            if (lhs[d] != rhs[d]) return false;
        return true;
    }
    friend constexpr bool operator!=(const Poly &lhs, const Poly &rhs) { return !(lhs == rhs); }

    friend decltype(auto) operator<<(std::ostream &lhs, const Poly &rhs) {
        auto s = 0, e = static_cast<int>(rhs.size());
        lhs << '[';
        for (auto &&i : rhs) {
            lhs << i;
            if (s >= 1) lhs << 'z';
            if (s > 1) lhs << '^' << s;
            if (++s != e) lhs << " + ";
        }
        if (e == 0) lhs << '0';
        return (lhs << ']');
    }
};

template<typename Iter>
Poly(Iter a, Iter b) -> Poly<typename std::iterator_traits<Iter>::value_type>;

}

#endif
