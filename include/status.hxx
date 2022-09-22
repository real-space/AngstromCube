#pragma once


  int constexpr STATUS_TEST_NOT_INCLUDED = -1;

#ifndef   STATUS_WITH_MESSAGE

  typedef int status_t;

#else  // STATUS_WITH_MESSAGE

#include <cstdio> // std::sprintf
#include <utility> // std::forward
#include <cstdint> // int32_t

class status_t
{
  private:
    int32_t _code;
    char _msg[124];
  public:

    status_t(int const code=0) : _code(code) {
        if (code) std::sprintf(_msg, "%i", code); else _msg[0] = 0;
    } // default constructor

    template <class... Args>
    status_t(char const *fmt, Args &&... args) : _code(1) {
        auto const nc = std::sprintf(_msg, fmt, std::forward<Args>(args)...);
        assert(nc <= 124); // out of bounds access
    } // constructor

    char const * message() const { return _code ? _msg : nullptr; }
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
  inline char const * status_message(status_t const & status) { return status.message(); }

#endif // STATUS_WITH_MESSAGE
