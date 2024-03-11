#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cassert> // assert
#include <cstdint> // uint32_t
#include <vector> // std::vector<T>
#include <algorithm> // std::fill

// #define debug_printf(...) std::printf(__VA_ARGS__)
#define debug_printf(...)

template <typename T>
class data_list { // a container for matrices with a variable number of columns per row
private:
    std::vector<T> data_; // raw data pointer to contiguous memory chunk [\sum_i m_[i]]
    std::vector<T*> ptrs_; // vector of pointers to starts [n_]
    std::vector<uint32_t> m_; // number of elements [n_]
    size_t mem_; // major memory consumption of data_
    uint32_t n_; // number of entries
    uint32_t max_m_; // largest entry of m_
public:

    template <typename int_t>
    data_list(uint32_t const n, int_t const ms[], T const init_value=T(0)) 
        : ptrs_(n, nullptr), m_(n), mem_(0), n_(n), max_m_(0) {
        assert(n == n_); // safety checks on upper limit of n
        size_t num{0}; // total number of elements
        for (uint32_t i = 0; i < n; ++i) {
            auto const m = uint32_t(std::max(ms[i], int_t(0)));
            m_[i] = m;
            assert(ms[i] == m_[i]); // safety checks on upper limit of ms[i]
            max_m_ = std::max(max_m_, m);
            num += m_[i];
        } // i
        mem_ = num*sizeof(T);
        debug_printf("# data_list() constructor tries to allocate %.6f MByte\n", mem_*1e-6);
        data_ = std::vector<T>(num, init_value);
        num = 0;
        for (uint32_t i = 0; i < n; ++i) {
            ptrs_[i] = &data_[num];
            num += m_[i];
        } // i
        assert(num*sizeof(T) == mem_); // consistency check
    } // constructor

    template <typename int_t>
    data_list(std::vector<int_t> const & ms, T const init_value=T(0)) 
      : data_list(ms.size(), ms.data(), init_value) {} // delegating constructor

    data_list(void) : data_(0), ptrs_(0), m_(0), mem_(0), n_(0), max_m_(0) {} // default constructor

    ~data_list() {
        debug_printf("# ~data_list() destructor tries to free %.6f MByte\n", mem_*1e-6);
    } // destructor

    data_list(data_list<T> && rhs) = delete; // move constructor
//  data_list(data_list<T> && rhs) {
//      debug_printf("# data_list(data_list<T> && rhs) { *this = std::move(rhs); }; \n");
//      *this = std::move(rhs);
//  } // move constructor

    data_list(data_list<T> const & rhs) = delete; // copy constructor
    data_list& operator= (data_list<T> const & rhs) = delete; // move assignment

    // move assignment is needed in many situations
//  data_list& operator= (data_list<T> && rhs) = delete; // move assignment
    data_list& operator= (data_list<T> && rhs) {
        debug_printf("# data_list& operator= (data_list<T> && rhs) transfers %.6f MByte\n", rhs.mem_*1e-6);
        data_.swap(rhs.data_);
        ptrs_.swap(rhs.ptrs_);
        m_.swap(rhs.m_);
        mem_   = rhs.mem_;    rhs.mem_   = 0;
        n_     = rhs.n_;      rhs.n_     = 0;
        max_m_ = rhs.max_m_;  rhs.max_m_ = 0;
        return *this;
    } // move assignment  


    // access operators
    T const & operator () (uint32_t const i, uint32_t const j) const { return ptrs_[i][j]; } // (i,j)
    T       & operator () (uint32_t const i, uint32_t const j)       { return ptrs_[i][j]; } // (i,j)

    T const & at(uint32_t const i, uint32_t const j) const { assert(i < n_); assert(j < m_[i]); return ptrs_[i][j]; }
    T       & at(uint32_t const i, uint32_t const j)       { assert(i < n_); assert(j < m_[i]); return ptrs_[i][j]; }

    T* operator[] (uint32_t const i) const { assert(i < n_); return ptrs_[i]; }

    // member access functions
    T const *const * data() const { return ptrs_.data(); } // read-only
    T       *const * data()       { return ptrs_.data(); }

    uint32_t nrows() const { return n_; } // number of rows
    uint32_t mcols() const { return max_m_; } // max. number of cols
    uint32_t ncols(uint32_t const i) const { assert(i < n_); return m_[i]; }
    uint32_t const * m() const { return m_.data(); } // numbers of cols
    size_t fill(T const value=T(0)) { std::fill(data_.begin(), data_.end(), value); return data_.size(); } // set value

}; // data_list

#undef debug_printf
