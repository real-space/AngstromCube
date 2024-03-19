#pragma once
// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, int16_t, int8_t
#include <cassert> // assert
#include <algorithm> // std::max
#include <vector> // std::vector<T>
#include <type_traits> // std::is_signed
#include <utility> // std::swap

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "print_tools.hxx" // printf_vector
#include "green_memory.hxx" // get_memory, free_memory
#include "simple_stats.hxx" // ::Stats<>

namespace green_sparse {

// #ifdef    DEBUG
    #define green_sparse_debug_printf(...) std::printf(__VA_ARGS__)
    int constexpr echo_copy = 9;
// #else  // DEBUG
//     #define green_sparse_debug_printf(...)
//     int constexpr echo_copy = 0;
// #endif // DEBUG

  char const int_t_names[][10] = {"uint8_t", "uint16_t", "uint32_t", "uint64_t", "uint128_t"};

  template <typename int_t>
  char const * int_t_name() {
      int const b = sizeof(int_t); // number of bytes
      int const log2b = (b < 2) ? 0 : ((b < 4) ? 1 : ((b < 8) ? 2 : ((b < 16) ? 3 : 4 ))); 
      return int_t_names[log2b] + std::is_signed<int_t>::value; // if signed, forward the char pointer by 1 to skip the 'u'
  } // int_t_name

  template <typename ColIndex_t=uint32_t, typename RowIndex_t=uint32_t>
  class sparse_t { // compressed sparse row (CSR) style sparsity descriptor
  public:

      static bool constexpr row_signed = std::is_signed<RowIndex_t>::value;
      static bool constexpr col_signed = std::is_signed<ColIndex_t>::value;

      sparse_t() : rowStart_(nullptr), colIndex_(nullptr), rowIndex_(nullptr), nRows_(0) {
          green_sparse_debug_printf("# sparse_t default constructor\n");
      } // default constructor

      sparse_t( // constructor
            std::vector<std::vector<ColIndex_t>> const & list
          , bool const with_rowIndex=false  // construct rowIndex now(true) or on demand(false)
          , char const *matrix_name=nullptr // name for log and debug output
          , int const echo=0 // log-level
      ) {
          auto const name = matrix_name ? matrix_name : "constructor";
          if (echo > 7) std::printf("# sparse_t<ColIndex_t=%s,RowIndex_t=%s> %s\n", int_t_name<ColIndex_t>(), int_t_name<RowIndex_t>(), name);
          nRows_ = list.size();
          assert(list.size() == nRows_ && "RowIndex_t cannot hold number of rows");
          auto const rowStart_nc = get_memory<RowIndex_t>(nRows_ + 1, echo, "rowStart_"); // _nc: non-const
          assert(rowStart_nc && "sparse_t failed to get_memory for rowStart");
          size_t nnz{0}; // number of non-zeros
          simple_stats::Stats<> st;
          rowStart_nc[0] = 0;
          for (RowIndex_t iRow = 0; iRow < nRows_; ++iRow) {
              auto const n = list[iRow].size(); // number of non-zero entires in this row
              st.add(n);
              rowStart_nc[iRow + 1] = rowStart_nc[iRow] + n; // create prefix sum
              nnz += n;
          } // iRow
          rowStart_ = rowStart_nc;
          assert(nnz == rowStart_[nRows_] && "Counting error");

          rowIndex_ = with_rowIndex ? rowIndex() : nullptr;

          auto const colIndex_nc = get_memory<ColIndex_t>(nnz, echo, "colIndex_"); // _nc: non-const
          assert(rowStart_ && "sparse_t failed to get_memory for colIndex");
          for (RowIndex_t iRow = 0; iRow < nRows_; ++iRow) {
              auto const n = list[iRow].size(); // number of non-zero entires in this row
              assert(rowStart_[iRow] + n == rowStart_[iRow + 1]);
              for (size_t j = 0; j < n; ++j) {
                  auto const jnz = rowStart_[iRow] + j;
                  colIndex_nc[jnz] = list[iRow][j];
              } // j
          } // iRow
          colIndex_ = colIndex_nc;

          if (echo > 7) std::printf("# sparse_t constructed with %d rows and %ld non-zero elements\n", uint32_t(nRows_), size_t(nnz));
          if (echo > 5) std::printf("# sparse_t %s columns per row stats [%g, %g +/- %g, %g]\n", name, st.min(), st.mean(), st.dev(), st.max());
          if (echo > 9) { 
              std::printf("# sparse_t.rowStart"); if (echo > 19) std::printf("(%p)", (void*)rowStart_); 
              std::printf("[0..%d]=", nRows_);  printf_vector(" %d", rowStart_, nRows_ + 1);
              std::printf("# sparse_t.colIndex"); if (echo > 19) std::printf("(%p)", (void*)colIndex_); 
              std::printf("[0..%ld-1]= ", nnz); printf_vector(" %d", colIndex_, nnz);
          } // echo
      } // constructor

      ~sparse_t() { // destructor
          green_sparse_debug_printf("# sparse_t destructor\n");
          if (rowStart_) free_memory(rowStart_);
          if (colIndex_) free_memory(colIndex_);
          if (rowIndex_) free_memory(rowIndex_);
          nRows_ = 0;
      } // destructor

   // sparse_t(sparse_t && rhs)                   = delete; // move constructor
      sparse_t(sparse_t && rhs) {                           // move constructor
          green_sparse_debug_printf("# sparse_t move constructor sparse_t<ColIndex_t=%s,RowIndex_t=%s>\n", int_t_name<ColIndex_t>(), int_t_name<RowIndex_t>());
//        green_sparse_debug_printf("# sparse_t.colIndex(%p)= ", (void*)rhs.colIndex()); printf_vector(" %d", rhs.colIndex(), rhs.nNonzeros());
          std::swap(this->rowStart_, rhs.rowStart_);
          std::swap(this->colIndex_, rhs.colIndex_);
          std::swap(this->rowIndex_, rhs.rowIndex_);
          std::swap(this->nRows_   , rhs.nRows_   );
//        green_sparse_debug_printf("# sparse_t.colIndex(%p)= ", (void*)colIndex()); printf_vector(" %d", colIndex(), nNonzeros());
      } // move constructor
      // always use std::move() to assign a sparse_t

 //   sparse_t & operator= (sparse_t && rhs)      = delete; // move assignment operator
      sparse_t & operator= (sparse_t && rhs) {              // move assignment operator
          green_sparse_debug_printf("# sparse_t move assignment sparse_t<ColIndex_t=%s,RowIndex_t=%s>\n", int_t_name<ColIndex_t>(), int_t_name<RowIndex_t>());
//        green_sparse_debug_printf("# sparse_t.colIndex(%p)= ", (void*)rhs.colIndex()); printf_vector(" %d", rhs.colIndex(), rhs.nNonzeros());
          std::swap(this->rowStart_, rhs.rowStart_);
          std::swap(this->colIndex_, rhs.colIndex_);
          std::swap(this->rowIndex_, rhs.rowIndex_);
          std::swap(this->nRows_   , rhs.nRows_   );
//        green_sparse_debug_printf("# sparse_t.colIndex(%p)= ", (void*)colIndex()); printf_vector(" %d", colIndex(), nNonzeros());
          return *this;
      } // move assignment operator

      sparse_t(sparse_t const & rhs)              = delete; // copy constructor
      sparse_t & operator= (sparse_t       & rhs) = delete; // copy assignment operator

 //   sparse_t & operator= (sparse_t const & rhs) = delete; // copy assignment operator
      sparse_t & operator= (sparse_t const & rhs) {         // copy assignment operator
          green_sparse_debug_printf("# sparse_t copy assignment sparse_t<ColIndex_t=%s,RowIndex_t=%s>\n", int_t_name<ColIndex_t>(), int_t_name<RowIndex_t>());
          int constexpr echo = echo_copy;
          nRows_ = rhs.nRows();
          auto const rowStart_nc = get_memory<RowIndex_t>(nRows_ + 1, echo, "rowStart_ (copy)");
          set(rowStart_nc, nRows_ + 1, rhs.rowStart()); // deep copy
          rowStart_ = rowStart_nc;
          auto const nnz = size_t(rhs.nNonzeros());
          auto const colIndex_nc = get_memory<ColIndex_t>(nnz, echo, "colIndex_ (copy)");
          set(colIndex_nc, nnz, rhs.colIndex()); // deep copy
          colIndex_ = colIndex_nc;
          return *this;
      } // copy assignment operator (deep copy)


#ifdef __NVCC__
      __host__ __device__ // the following member function is also called on GPUs
#endif
      RowIndex_t const * rowStart() const { return rowStart_; };
#ifdef __NVCC__
      __host__ __device__ // the following member function is also called on GPUs
#endif
      ColIndex_t const * colIndex() const { return colIndex_; };
      RowIndex_t            nRows() const { return nRows_; }
      RowIndex_t        nNonzeros() const { return (row_signed && nRows_ < 0) ? 0 : (rowStart_ ? rowStart_[nRows_] : 0); }

  private:

      // __host__ 
      bool invalid_row_index_(RowIndex_t const iRow) const {
          if (iRow >= nRows_)  return true;  // invalid
          else if (iRow < 0)   return true;  // invalid
          else                 return false; //   valid
      } // invalid_row_index_

  public:

#if 0 // are these methods ever needed?

      RowIndex_t nNonzeroCols(RowIndex_t const iRow) const {
          if (invalid_row_index_(iRow)) return 0;
          return rowStart_ ? rowStart_[iRow] : 0;
      } // nNonzeroCols

      ColIndex_t const * nonzeroCols(RowIndex_t const iRow) const {
          if (invalid_row_index_(iRow)) return nullptr;
          auto const row_start = rowStart_ ? rowStart_[iRow] : 0;
          return colIndex_ ? &colIndex_[row_start] : nullptr;
      } // nonzeroCols

#endif // 0

      bool is_in(RowIndex_t const iRow, ColIndex_t const jCol, RowIndex_t *index=nullptr) const {
          if (invalid_row_index_(iRow)) return false;
          if (nullptr == rowStart_ || nullptr == colIndex_) return false;
//        std::printf("# search for index (iRow=%d, jCol=%d)\n", iRow, jCol);
          for (auto jnz = rowStart_[iRow]; jnz < rowStart_[iRow]; ++jnz) {
              if (jCol == colIndex_[jnz]) {
                  if (index) *index = jnz; // export non-zero index (optional)
                  return true;
              } // match
          } // jnz
          return false;
      } // is_in

      // needed?
      RowIndex_t const * rowIndex() { // non-const member function triggers the generation of internal rowIndex_
          if (nullptr == rowIndex_) {
              // try to construct the rowIndex list
              auto const nnz = nNonzeros();
              if (nnz < 1) return nullptr;
              auto const rowIndex_nc = get_memory<RowIndex_t>(nnz); // _nc: non-const
              for (RowIndex_t iRow = 0; iRow < nRows_; ++iRow) {
                  for (auto jnz = rowStart_[iRow]; jnz < rowStart_[iRow + 1]; ++jnz) {
                      rowIndex_nc[jnz] = iRow;
                  } // jnz
              } // iRow
              rowIndex_ = rowIndex_nc;
          } // nullptr
          return rowIndex_;
      } // rowIndex

  private: // members
      RowIndex_t const *rowStart_ = nullptr;
      ColIndex_t const *colIndex_ = nullptr;
      RowIndex_t const *rowIndex_ = nullptr; // optional, will be constructed on first demand
      RowIndex_t nRows_ = 0;

  }; // sparse_t

  status_t all_tests(int const echo=0); // declaration only

#undef green_sparse_debug_printf

} // namespace green_sparse

