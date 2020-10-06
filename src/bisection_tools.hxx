#pragma once

#include <cstdio> // printf
#include <cassert> // assert

#include "status.hxx" // status_t

namespace bisection_tools {

  template<typename real_t=double>
  class bisector_t
  {
      private:
        real_t _x[2];
        real_t _y[2];
        real_t _thres;
        real_t _sign;
        int  _iter;
        int  _maxiter;
        char _mode;
        bool _is_converged;

        real_t suggest_x() const { return ('*' == _mode) ? std::sqrt(_x[0] * _x[1]) : (_x[0] + _x[1])/2; }
        
      public:
        
        bisector_t(real_t const range_begin, real_t const range_end
             , real_t const threshold=1e-6, char const mode='+', int const sign=1, int const maxiter=99)
          : _thres(threshold), _sign(sign), _iter(-2), _maxiter(maxiter), _mode(mode), _is_converged(false) {
            _x[0] = range_begin; _x[1] = range_end;
        } // constructor

        
        bool root(real_t & new_x, real_t const last_y, int const echo=0) {
            if (_is_converged) return false;
            if (_iter > _maxiter) return false;

            auto const last_x = new_x; // get a copy before new_x is modified
            if (-2 == _iter) {
                // ignore last_y since that cannot contain any value before the first iteration
                new_x = _x[0];
            } else if (-1 == _iter) {
                _y[0] = last_y;
                if (echo > 3) printf("# %s range_begin x[0]= %g\ty[0]= %g\n", __func__, last_x, last_y);    
                new_x = _x[1];
            } else {
                if (0 == _iter) {
                    if (echo > 3) printf("# %s range_end   x[1]= %g\ty[1]= %g\n", __func__, last_x, last_y);      
                    _y[1] = last_y;
                    if (_y[1]*_y[0] > 0) {
                        if (echo > 2) printf("# %s in zero-search y-values need to have opposite sign\n", __func__);      
                        return false; // failed
                    }
                    // determine the sub-type: 1:ascending -1:descending
                    if ((_sign*_y[1] < 0) && (_sign*_y[0] > 0)) {
                        _sign = -_sign; // is descending
                        if (echo > 3) printf("# %s toggle between ascending/descending\n", __func__);      
                    }
                } else {
                    if (echo > 7) printf("# %s iteration #%i  x= %g\ty= %g\n", __func__, _iter, last_x, last_y);
                } // iteration zero

                int const i01 = (_sign*last_y > 0); // i01 is the index which to overwrite
                _x[i01] = last_x;
                _y[i01] = last_y;

                _is_converged = (std::abs(_x[1] - _x[0]) < _thres); // returns false in the next iteration

                // suggest a new x-value
                new_x = suggest_x();
//              printf("# %s suggest x= %g\n", __func__, new_x);                
            }
            ++_iter;
            return true; // run again
        } // root

        bool extremum(real_t & new_x, real_t const last_y, int const echo=0) {
            if (_is_converged) return false;
            if (_iter > _maxiter) return false;

            auto const last_x = new_x; // get a copy before new_x is modified
            if (-2 == _iter) {
                // ignore last_y since that cannot contain any value before the first iteration
                new_x = _x[0];
            } else if (-1 == _iter) {
                _y[0] = last_y;
                if (echo > 3) printf("# %s range_begin x[0]= %g\ty[0]= %g\n", __func__, last_x, last_y);    
                new_x = _x[1];
            } else if (0 == _iter) {
                _y[1] = last_y;
                if (echo > 3) printf("# %s range_end   x[1]= %g\ty[1]= %g\n", __func__, last_x, last_y);
                new_x = suggest_x();
            } else {
                if (echo > 7) printf("# %s iteration #%i  x= %g\ty= %g\n", __func__, _iter, last_x, last_y);
                if (echo > 7) printf("# %s iteration #%i  y[0]= %g  y= %g  y[1]= %g\n", __func__, _iter, _y[0], last_y, _y[1]);
                int const i01 = (_sign*_y[0] > _sign*_y[1]); // i01 is the index which to overwrite
                if (echo > 7) printf("# %s iteration #%i  overwrite x[%i]= %g while x[%i]= %g\n", __func__, _iter, i01, _x[i01], 1-i01, _x[1-i01]);
                _x[i01] = last_x;
                _y[i01] = last_y;

                _is_converged = (std::abs(_x[1] - _x[0]) < _thres); // returns false in the next iteration

                // suggest a new x-value
                new_x = suggest_x();
//              printf("# %s suggest x= %g\n", __func__, new_x);                
            }
            ++_iter;
            return true; // run again
        } // extremum
        
        
        int get_iterations_needed() const { return _is_converged ? _iter : -_iter; }
        int is_converged() const { return _is_converged; }

  }; // bisector_t
  
  
#ifdef  NO_UNIT_TESTS
  inline status_t all_tests(int const echo=0) { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS
  
  inline status_t test_bisection_zero(int const echo=3) {
      bisector_t<double> b(0., 3., 1e-15);
      double x, y{0};
      while(b.root(x, y, echo)) {
          y = std::cos(x);
      } // while bisection
      if (echo > 0) printf("# %s solution y = sin(x= %.15f) = %g needed %d iterations\n",
                              __func__, x, y, b.get_iterations_needed());
      return 0;
  } // test_bisection_zero

  inline status_t test_bisection_maximum(int const echo=3) {
      bisector_t<double> b(0., 3., 1e-15);
      double x, y{0};
      while(b.extremum(x, y, echo)) {
          y = std::sin(x);
      } // while bisection
      if (echo > 0) printf("# %s solution y = sin(x= %.15f) = %g needed %d iterations\n",
                              __func__, x, y, b.get_iterations_needed());
      return 0;
  } // test_bisection_maximum

  inline status_t all_tests(int const echo=0) {
      if (echo > 0) printf("\n# %s %s\n", __FILE__, __func__);
      status_t status(0);
      status += test_bisection_zero(echo);
      status += test_bisection_maximum(echo);
      return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace bisection_tools
