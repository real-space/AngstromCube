#include <cstdint> // uint64_t
#include <string> // std::string
#include <cstring> // std::strrchr
#include <cstdio> // std::printf, ::snprintf
#include <cassert> // assert
#include <map> // std::map
#include <utility> // std::pair<T1, T2>, ::make_pair

#include "recorded_warnings.hxx" // MaxMessageLength

namespace recorded_warnings {

  // a simple hash function maps only 9 out of 216000 english words onto the same integer number
  inline uint64_t simple_string_hash(char const *string) {
      uint64_t hash{5381};
      char const *c = string;
      while (*c != 0) {
         hash = hash * 33 + uint64_t(*c);
         ++c;
      } // while
      return hash;
  } // simple_string_hash

  inline uint64_t combined_hash(char const *file, int const line) {
      int const LineBits = 14; // host the line number within the first 14 bits
      return (simple_string_hash(file) << LineBits) | (line & ((1ul << LineBits) - 1));
  } // combinded_hash

  class WarningRecord {
  private:
      char*       message_;
      uint64_t    hash_;
      std::string source_file_name_;
      std::string function_name_;
      uint64_t    times_overwritten_;
      uint32_t    times_printed_;
      uint32_t    source_file_line_;
      // ToDo: we can include the time of last overwrite, but then, log files differ between two executions
  public:

      WarningRecord(char const *file, int const line, char const *func=nullptr)
        : message_(new char[MaxMessageLength]) // this memory is never released again
        , hash_(combined_hash(file, line))
        , source_file_name_(file)
        , function_name_(func)
        , times_overwritten_(0)   
        , times_printed_(0)
        , source_file_line_(line)
      {
#ifdef  DEBUG
          std::printf("# WarningRecord:constructor allocates a new warning message string with"
              " max. %d chars at %p\n# ... for warnings launched at %s:%d --> hash = %16llx\n",
                     message_length, (void*)message_,             file,line,  hash_);
#endif // DEBUG
      } // constructor

      ~WarningRecord(void) {
#ifdef  DEBUG          
          std::printf("# WarningRecord:destructor: old warning message"
                      " at %p for warnings launched at %s:%d reads:\n#\t%s\n", 
                      (void*)message_, get_sourcefile(), source_file_line_, message_);
#endif // DEBUG
          // delete[] message_; // seems like this happens automagically
      } // destructor

      char* get_message(void) { ++times_overwritten_; return message_; }
      char* get_message_pointer(void) const { return message_; }
      char const* get_sourcefile(void) const { return source_file_name_.c_str(); }
      char const* get_functionname(void) const { return function_name_.c_str(); }
      int get_sourceline(void) const { return source_file_line_; }
      size_t get_times(void) const { return times_overwritten_; }
      int get_times_printed(void) const { return times_printed_; }
      void increment_times_printed(void) { ++times_printed_; }

  }; // class WarningRecord



  std::pair<char*,int> _manage_warnings(char const *file, int const line, char const *func, int const echo=0) {
    if (echo > 8) std::printf("\n# %s:%d  %s(file=%s, line=%d, echo=%d)\n", 
                      __FILE__, __LINE__, __func__, file, line, echo);

    static std::map<uint64_t, WarningRecord> map_;
    if (line < 1) { // line numbers created by the preprocessor start from 1
        assert('?' == file[0]); // make sure that we want special functionality
        
        if (0 == line) {
            // show_warnings() has been called
            if (echo > 0) {
                auto const nw = map_.size();
                if ((echo < 3) || (nw < 1)) { 
                    // only give a summary of how many
                    std::printf("# %ld warnings have been recorded.\n", nw);
                } else {
                    std::printf("\n#\n# recorded %ld different warnings:\n", nw);
                    size_t total_count{0};
                    for (auto &hw : map_) {
                        auto const &w = hw.second;
                        auto const n_times = w.get_times();
                        std::printf("# \tin %s:%d %s (%ld times)\n"
                               "# \t\t%s\n", w.get_sourcefile(), 
                            w.get_sourceline(), w.get_functionname(), 
                            n_times, w.get_message_pointer());
                        total_count += n_times;
                    } // w
                    if (nw > 0) std::printf("# %ld warnings in total\n", total_count);
                } // summary
            } // echo
        } else { // line
            // clear_warnings() has been called
            if (echo > 1) std::printf("# clear all %ld warnings from records\n", map_.size());
            map_.clear();
        } // line
        return std::make_pair(nullptr, 0);

    } else { // special functions
      
        // regular usage with file and line
      
        // remove path name to source file
        auto const short_file = after_last_slash(file);
        auto const hash = combined_hash(short_file, line);
        
        WarningRecord *w;
        auto const search = map_.find(hash);
        if (map_.end() != search) {
            if (echo > 1) std::printf("# %s: found entry for hash %16llx\n", __func__, hash);
            w = &search->second;
        } else {
            if (echo > 1) std::printf("# %s: insert new entry for hash %16llx\n", __func__, hash);
            auto const iit = map_.insert({hash, WarningRecord(short_file, line, func)});
            w = &iit.first->second;
        } // found

        // output the warning to stdout and stderr when encountered the 1st time, otherwise,
        // we could have a segfault later and do not know where that could be coming from
        
        // configuration:
        int const MaxWarningsToErr =  1; // show max M warnings in stderr, negative means all
        int const MaxWarningsToLog =  3; // show max |M| warnings in stdout, negative
                // values suppress the "# This warning will not be shown again!" info

        int flags{0};
        int const times_printed = w->get_times_printed();
        if (times_printed < std::abs(MaxWarningsToLog)) {
            flags |= 1; // 1: message to stdout
            if (times_printed == MaxWarningsToLog - 1) 
                flags |= 4; // "# This warning will not be shown again!"
        } // print to stdout
        if ((MaxWarningsToErr < 0) || (times_printed < MaxWarningsToErr)) flags |= 2; // 2: message to stderr
        if (flags) w->increment_times_printed(); // will print message to stdout or stderr

        return std::make_pair(w->get_message(), flags);
    } // special functions

  } // _manage_warnings

  std::pair<char*,int> _new_warning(char const *file, int const line, char const *func=nullptr) {
      assert(line >= 0 && "line numbers created by the preprocessor start from 1");
      return _manage_warnings(file, line, func);
  } // _new_warning

  status_t show_warnings(int const echo) {
      return (nullptr != _manage_warnings("?",  0, nullptr, echo).first);
  } // show_warnings

  status_t clear_warnings(int const echo) {
      return (nullptr != _manage_warnings("?", -1, nullptr, echo).first);
  } // clear_warnings


#ifdef  NO_UNIT_TESTS
  status_t all_tests(int const echo) { return STATUS_TEST_NOT_INCLUDED; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
      if (echo > 1) std::printf("\n# %s:%d  %s\n\n", __FILE__, __LINE__, __func__);
      WarningRecord wr(__FILE__,__LINE__,__func__);
      auto const msg = wr.get_message();
      auto const nchars = std::snprintf(msg, MaxMessageLength, 
            "This is a non-recorded warning from %s:%d", __FILE__, __LINE__);
      return (nchars >= MaxMessageLength);
  } // test_create_and_destroy

  status_t test_preprocessor_macro(int const echo=9) {
      if (echo > 1) std::printf("\n# %s:%d  %s\n", __FILE__, __LINE__, __func__);
      auto const nchars = warn("This is a test warning from %s:%d", __FILE__, __LINE__);
      return (nchars < 30); // error if the warning is not printed
  } // test_preprocessor_macro

  status_t test_overwriting(int const echo=9) {
      if (echo > 1) std::printf("\n# %s:%d  %s\n", __FILE__, __LINE__, __func__);
      int nchars{0};
      for (int i = 0; i < 9; ++i) {
          nchars = warn("This is a test warning from inside a loop, iteration #%d", i);
          if (echo > 19) std::printf("# %s: warning message had %d characters\n", __func__, nchars);
      } // i
      return nchars; // returns 0 if the warnings were muted before the 9th iteration
  } // test_overwriting

  status_t all_tests(int const echo) {
      status_t stat(0);
      stat += test_create_and_destroy(echo);
      stat += test_preprocessor_macro(echo);
      stat += test_overwriting(echo);
      // clean up
      stat += show_warnings(echo); // display those warnings launched for test purposes
      stat += clear_warnings(echo); // clear test warnings from record, this is necessary
                                    // so they do not appear when running all unit tests
      return stat;
  } // all_tests

#endif // NO_UNIT_TESTS  

} // namespace recorded_warnings
