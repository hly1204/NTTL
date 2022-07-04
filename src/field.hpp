#ifndef FIELD_HPP
#define FIELD_HPP

#include <concepts>
#include <cstdint>
#include <utility>
#include <type_traits>

namespace nttl {

namespace detail {

constexpr auto is_prime(std::uint32_t v) {
    if (v <= 1) return false;
    for (auto i = 2u; static_cast<std::uint64_t>(i) * i <= v; ++i)
        if (v % i == 0) return false;
    return true;
}

}

template<std::uint32_t V>
static inline constexpr auto is_prime_v = detail::is_prime(V);

template<typename T>
concept Field =
    std::default_initializable<T> &&
    requires(T u, const T &v) {
        {v.inv()} -> std::convertible_to<T>;
        {-v} -> std::same_as<T>;
        {u += v} -> std::same_as<std::add_lvalue_reference_t<T>>;
        {u -= v} -> std::same_as<std::add_lvalue_reference_t<T>>;
        {u *= v} -> std::same_as<std::add_lvalue_reference_t<T>>;
        {u /= v} -> std::convertible_to<std::add_lvalue_reference_t<T>>;
        {u + v} -> std::same_as<T>;
        {u - v} -> std::same_as<T>;
        {u * v} -> std::same_as<T>;
        {u / v} -> std::convertible_to<T>;
        {u ^ -1ll} -> std::convertible_to<T>; // 同 pow
        {u.pow(-1ll)} -> std::convertible_to<T>;
        {u == T{}} -> std::same_as<bool>;
    };

template<typename T>
concept FiniteField =
    Field<T> &&
    std::is_integral_v<std::decay_t<decltype(T::x())>> &&
    (std::is_same_v<decltype(T::card().first), std::uint32_t> ||
    std::is_same_v<decltype(T::card().first), std::uint64_t>) &&
    std::is_same_v<decltype(T::card().second), std::uint32_t> &&
    requires(T v) {
        v[0];
        requires Field<std::decay_t<decltype(v[0])>>;
    } &&
    requires(const T v) {
        v[0];
        requires Field<std::decay_t<decltype(v[0])>>;
    };
}

#endif
