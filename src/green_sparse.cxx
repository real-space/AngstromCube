// This file is part of AngstromCube under MIT License

#include <cstdio> // std::printf
#include <cstdint> // int64_t, int32_t, int16_t, int8_t
#include <vector> // std::vector<T>

#define DEBUG

#include "green_sparse.hxx"

#include "status.hxx" // status_t, STATUS_TEST_NOT_INCLUDED
#include "print_tools.hxx" // printf_vector

namespace green_sparse {

#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo=0) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_sparse(int const echo=0, int const n=7) {
      sparse_t<> s0;              // test default constructor
      sparse_t<int> sint;         // test default constructor
      sparse_t<int,int> sintint;  // test default constructor
      std::vector<std::vector<uint8_t>> list(n);
      for (int i = 0; i < n; ++i) {
          list[i].resize(i + 1); // lower triangular matrix
          for (int j = 0; j <= i; ++j) list[i][j] = j;
      } // i
      sparse_t<uint8_t> s(list, false, "lowerTriangularMatrix", echo); // test constructor
      if (echo > 3) std::printf("# %s nRows= %d, nNonzeros= %d\n", __func__, s.nRows(), s.nNonzeros());
      if (echo > 6) std::printf("# %s last element = %d\n", __func__, int(s.colIndex()[s.nNonzeros() - 1]));
      if (echo > 6) { std::printf("# rowIndex "); printf_vector(" %d", s.rowIndex(), s.nNonzeros()); }
      if (echo > 3) std::printf("# sizeof(sparse_t) = %ld Byte\n", sizeof(sparse_t<>));
      return 0;
  } // test_sparse

  status_t test_experiment(int const echo=0, int const n=3, int const m=4) {
      if (echo > 3) std::printf("\n# %s:\n", __func__);
      // sparse_t<> s0; // test default constructor
      // auto s1(s0); // copy constructor is deleted
      std::vector<std::vector<uint32_t>> list(n, std::vector<uint32_t>(m, 3));
      sparse_t<> s2(list, false, "test_s2", echo);
      auto s3 = sparse_t<>(list, false, "test_s3", echo);
      return 0;
  } // test_experiment

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_sparse(echo);
      stat += test_experiment(echo);
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace green_sparse
