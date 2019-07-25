#pragma once

// template<int nbits> size_t inline align(int64_t const in) 
//   { return ((((in - 1) >> nbits) + 1) << nbits); }
template<int nbits> size_t inline align(int64_t const in) { return in; }
