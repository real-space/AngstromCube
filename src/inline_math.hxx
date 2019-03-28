#pragma once

template <typename T> inline int sgn(T const val) {
    return (T(0) < val) - (val < T(0));
} // sgn

template <typename T> inline T pow2(T const x) { return x*x; } // pow2
