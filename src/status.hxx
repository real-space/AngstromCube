#pragma once

// #define STATUS_WITH_MESSAGE

#ifndef STATUS_WITH_MESSAGE
  typedef int status_t;
  
  template <class... Args>
  set_status(char const *fmt, Args &&... args) { return 1; }
  set_status(int const i) { return i; }
#else

#include <cstdio> // printf, std::sprintf
#include <utility> // std::forward
#include <cstdint> // int32_t

class status_t 
{
  private:
    int32_t _code;
    char _msg[124];
  public:

    status_t(int const code=0) : _code(code) { 
        _msg[0] = 0; // empty message
    } // default constructor
    
    template <class... Args>
    status_t(char const *fmt, Args &&... args) : _code(1) {
        std::sprintf(_msg, fmt, std::forward<Args>(args)...);
    } // constructor

    char const * message() const { return _msg; }
    int  const      code() const { return _code; }
    
    // we want to ask if(status)
    operator bool() const { return (0 != _code); };

    // we want to be able to add status variables
    operator int() const { return _code; }
    
    status_t & operator += (int const & rhs) { _code += rhs; return *this; }
    status_t & operator ++() { ++_code; return *this; }
    bool operator > (int const rhs) { return _code > rhs; }
    
}; // class

#define set_status status_t

#endif // STATUS_WITH_MESSAGE
