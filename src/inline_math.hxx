#pragma once

template <typename T> inline int sgn(T const val) {
    return (T(0) < val) - (val < T(0));
} // sgn

template <typename T> inline T pow2(T const x) { return x*x; } // pow2
template <typename T> inline T pow3(T const x) { return x*pow2(x); } // pow3

template <typename real_t> 
inline real_t intpow(real_t const x, unsigned n) {
	real_t xbin = x, xpow = 1;
	while (n) {
		if (n & 1) xpow *= xbin; // if n modulo 2 == 1
		n >>= 1; // divide n by 2
		xbin *= xbin; // square x
	}
	return xpow;
} // intpow
