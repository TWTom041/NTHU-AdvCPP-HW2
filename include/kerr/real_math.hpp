// real_math.hpp - uniform math wrappers so that the Kerr engine can be
// instantiated with float / double / long double / __float128 / std::float16_t.
//
// Every function the engine needs lives in namespace kerr::m and is overloaded
// per type, so templates simply call m::sin(x) etc. (std:: does not know
// __float128, and ADL does not work for builtin types.)
#pragma once

#include <cmath>
#include <limits>
#include <string>
#ifdef __SIZEOF_FLOAT128__  // x86-64: binary128 via libquadmath (aarch64: long double is binary128)
#define KERR_HAS_FLOAT128 1
#include <quadmath.h>
#endif
#if __has_include(<experimental/simd>)
#include <experimental/simd>
#endif
#if __has_include(<stdfloat>)
#include <stdfloat>
#endif

namespace kerr::m {

#define KERR_STD_OVERLOADS(T)                                              \
    inline T sqrt(T x) { return std::sqrt(x); }                            \
    inline T sin(T x) { return std::sin(x); }                              \
    inline T cos(T x) { return std::cos(x); }                              \
    inline T acos(T x) { return std::acos(x); }                            \
    inline T log(T x) { return std::log(x); }                              \
    inline T fabs(T x) { return std::fabs(x); }                            \
    inline bool isfinite(T x) { return std::isfinite(x); }

KERR_STD_OVERLOADS(float)
KERR_STD_OVERLOADS(double)
KERR_STD_OVERLOADS(long double)
#undef KERR_STD_OVERLOADS

#ifdef KERR_HAS_FLOAT128
// IEEE binary128, implemented in software by libgcc (__addtf3, __multf3 ...)
// and libquadmath (sinq, cosq ...).
inline __float128 sqrt(__float128 x) { return sqrtq(x); }
inline __float128 sin(__float128 x) { return sinq(x); }
inline __float128 cos(__float128 x) { return cosq(x); }
inline __float128 acos(__float128 x) { return acosq(x); }
inline __float128 log(__float128 x) { return logq(x); }
inline __float128 fabs(__float128 x) { return fabsq(x); }
inline bool isfinite(__float128 x) { return finiteq(x); }
#endif

#ifdef __STDCPP_FLOAT16_T__
// IEEE binary16: math functions are evaluated in float and rounded back.
using f16 = std::float16_t;
inline f16 sqrt(f16 x) { return f16(std::sqrt(float(x))); }
inline f16 sin(f16 x) { return f16(std::sin(float(x))); }
inline f16 cos(f16 x) { return f16(std::cos(float(x))); }
inline f16 acos(f16 x) { return f16(std::acos(float(x))); }
inline f16 log(f16 x) { return f16(std::log(float(x))); }
inline f16 fabs(f16 x) { return f16(std::fabs(float(x))); }
inline bool isfinite(f16 x) { return std::isfinite(float(x)); }
#endif

#if __has_include(<experimental/simd>)
// SIMD packets (std::experimental::simd, Parallelism TS v2) so that
// Spacetime<simd>::rhs can advance several rays in lock-step.
namespace stdx = std::experimental;
template <class T, class A> inline stdx::simd<T, A> sqrt(const stdx::simd<T, A>& x) { return stdx::sqrt(x); }
template <class T, class A> inline stdx::simd<T, A> sin(const stdx::simd<T, A>& x) { return stdx::sin(x); }
template <class T, class A> inline stdx::simd<T, A> cos(const stdx::simd<T, A>& x) { return stdx::cos(x); }
template <class T, class A> inline stdx::simd<T, A> fabs(const stdx::simd<T, A>& x) { return stdx::abs(x); }
#endif

#ifdef KERR_HAS_FLOAT128
template <class T> constexpr T pi() { return T(3.14159265358979323846264338327950288419716939937510582Q); }
#else
template <class T> constexpr T pi() { return T(3.14159265358979323846264338327950288419716939937510582L); }
#endif

template <class T> inline T min(T a, T b) { return b < a ? b : a; }
template <class T> inline T max(T a, T b) { return a < b ? b : a; }

}  // namespace kerr::m

namespace kerr {

template <class T> const char* type_name();
template <> inline const char* type_name<float>() { return "float"; }
template <> inline const char* type_name<double>() { return "double"; }
template <> inline const char* type_name<long double>() { return "long double"; }
#ifdef KERR_HAS_FLOAT128
template <> inline const char* type_name<__float128>() { return "__float128"; }
#endif
#ifdef __STDCPP_FLOAT16_T__
template <> inline const char* type_name<std::float16_t>() { return "float16_t"; }
#endif

// Machine epsilon for every supported type (numeric_limits<__float128> is
// only specialised in GNU mode, so do it by hand).
template <class T> inline T epsilon() { return std::numeric_limits<T>::epsilon(); }
#ifdef KERR_HAS_FLOAT128
template <> inline __float128 epsilon<__float128>() { return FLT128_EPSILON; }
#endif

}  // namespace kerr
