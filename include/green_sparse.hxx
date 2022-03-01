#pragma once

#include <cstdint> // int64_t, int32_t, int16_t, int8_t
#include <cassert> // assert
#include <algorithm> // std::max
#include <utility> // std::swap
#include <vector> // std::vector<T>
#include <cstdio> // std::printf
#include <type_traits> // std::is_signed

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "print_tools.hxx" // printf_vector
#include "green_memory.hxx" // get_memory, free_memory
#include "simple_stats.hxx" // ::Stats<>

#define DEBUG

#ifdef  DEBUG
  #define debug_printf(...) std::printf(__VA_ARGS__)
#else  // DEBUG
  #define debug_printf(...)
#endif // DEBUG

namespace green_sparse {

  double const MByte = 1e-6; char const *const _MByte = "MByte";

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

      sparse_t() : rowStart_(nullptr), colIndex_(nullptr), rowIndex_(nullptr), nRows_(0) {
          debug_printf("# sparse_t default constructor\n");
      } // default constructor

      sparse_t(
            std::vector<std::vector<ColIndex_t>> const & list
          , bool const with_rowIndex=false  // construct rowIndex now(true) or on demand(false)
          , char const *matrix_name=nullptr // name for log and debug output
          , int const echo=0 // log-level
      ) {
          char const *name = matrix_name ? matrix_name : "constructor";
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
          if (echo > 9) { std::printf("# sparse_t.rowStart(%p)= ", (void*)rowStart_); printf_vector(" %d", rowStart_, nRows_); }
          if (echo > 9) { std::printf("# sparse_t.colIndex(%p)= ", (void*)colIndex_); printf_vector(" %d", colIndex_, nnz); }
      } // constructor

      ~sparse_t() {
          debug_printf("# sparse_t destructor\n");
          if (rowStart_) free_memory(rowStart_);
          if (colIndex_) free_memory(colIndex_);
          if (rowIndex_) free_memory(rowIndex_);
          nRows_ = 0;
      } // destructor

      sparse_t(sparse_t && rhs) = delete;

      sparse_t(sparse_t const & rhs) = delete;

      sparse_t & operator= (sparse_t       & rhs) = delete; // copy assignment operator
      sparse_t & operator= (sparse_t const & rhs) = delete; // copy assignment operator

      sparse_t & operator= (sparse_t && rhs) {
          debug_printf("# sparse_t move assignment sparse_t<ColIndex_t=%s,RowIndex_t=%s>\n", int_t_name<ColIndex_t>(), int_t_name<RowIndex_t>());
//        debug_printf("# sparse_t.colIndex(%p)= ", (void*)rhs.colIndex()); printf_vector(" %d", rhs.colIndex(), rhs.nNonzeros());
          std::swap(rowStart_, rhs.rowStart_);
          std::swap(colIndex_, rhs.colIndex_);
          std::swap(rowIndex_, rhs.rowIndex_);
          std::swap(nRows_   , rhs.nRows_);
//        debug_printf("# sparse_t.colIndex(%p)= ", (void*)colIndex()); printf_vector(" %d", colIndex(), nNonzeros());
          return *this;
      } // move assignment

      __host__ __device__ RowIndex_t const * rowStart() const { return rowStart_; };
      __host__ __device__ ColIndex_t const * colIndex() const { return colIndex_; };
      __host__ __device__ RowIndex_t            nRows() const { return nRows_; }
      __host__ __device__ RowIndex_t        nNonzeros() const { return (nRows_ < 0) ? 0 : (rowStart_ ? rowStart_[nRows_] : 0); }

  private:
    
      __host__ bool invalid_row_index_(RowIndex_t const iRow) const {
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

      __host__ bool is_in(RowIndex_t const iRow, ColIndex_t const jCol, RowIndex_t *index=nullptr) const {
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

      __host__ RowIndex_t const * rowIndex() { // non-const member function triggers the generation of internal rowIndex_
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
      RowIndex_t const *rowStart_;
      ColIndex_t const *colIndex_;
      RowIndex_t const *rowIndex_; // optional, will be constructed on first demand
      RowIndex_t nRows_;

  }; // sparse_t



#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  inline status_t test_sparse(int const echo=0, int const n=7) {
      sparse_t<> s0; // test default constructor
      sparse_t<int> sint; // test default constructor
      sparse_t<int,int> sintint; // test default constructor
      std::vector<std::vector<uint8_t>> list(n);
      for (int i = 0; i < n; ++i) {
          list[i].resize(i + 1); // lower triangular matrix
          for (int j = 0; j <= i; ++j) list[i][j] = j;
      } // i
      sparse_t<uint8_t> s(list, false, "lowerTriangularMatrix"); // test constructor
      if (echo > 3) std::printf("# %s nRows= %d, nNonzeros= %d\n", __func__, s.nRows(), s.nNonzeros());
      if (echo > 6) std::printf("# %s last element = %d\n", __func__, int(s.colIndex()[s.nNonzeros() - 1]));
      if (echo > 6) { std::printf("# rowIndex "); printf_vector(" %d", s.rowIndex(), s.nNonzeros()); }
      if (echo > 3) std::printf("# sizeof(sparse_t) = %ld Byte\n", sizeof(sparse_t<>));
      return 0;
  } // test_sparse

  inline status_t all_tests(int const echo=0) {
      status_t stat(0);
      stat += test_sparse(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_sparse

#undef debug_printf
