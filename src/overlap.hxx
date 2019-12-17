typedef int status_t;

namespace overlap {

  template<typename real_t>
  status_t generate_tensor(real_t tensor[], int const ncut, 
                     double const sigma=2, // 2:typical for density tensor
                     double const sigma0=1, double const sigma1=1);
  
  status_t all_tests();

} // namespace overlap
