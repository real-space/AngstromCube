typedef int status_t;

namespace angular_grid {

  int Lebedev_grid_size(int const ellmax, int echo=0);
  
  template<typename real_t>
  status_t create_Lebedev_grid(int const ellmax, real_t xyzw[][4], int echo=9);

  status_t all_tests();
  
} // namespace angular_grid
