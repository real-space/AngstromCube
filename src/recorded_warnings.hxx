#pragma once

#include <cstdint> // uint64_t
#include <string> // std::string
#include <cstdio> // printf

typedef int status_t;

#define warn recorded_warnings::new_warning(__FILE__, __LINE__) 

namespace recorded_warnings {

  char* new_warning(char const *file, int const line);

  status_t show_warnings(int const echo=1);

  status_t clear_warnings(int const echo=1);

  status_t all_tests();

} // namespace recorded_warnings
