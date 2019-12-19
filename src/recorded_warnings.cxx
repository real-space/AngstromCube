#include <cstdint> // uint64_t
#include <string> // std::string
#include <cstdio> // printf, sprintf
#include <cassert> // assert
#include <map> // std::map

#include "recorded_warnings.hxx"

namespace recorded_warnings {

  // a simple hash function maps only 9 out of 216000 english words onto the same integer number
  inline uint64_t simple_string_hash(char const *string) {
      uint64_t hash = 5381;
      char const *c = string;
      while (*c != 0) {
         hash = hash * 33 + (*c);
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
      char*     message_;
      uint64_t  hash_;
      std::string source_file_name_;
      uint64_t    times_overwritten_;
      uint32_t    source_file_line_;
      // ToDo: we can include the time of last overwrite
  public:
      static const int DefaultMessageLength = 400;

      WarningRecord(char const *file, int const line, 
                  int const message_length=DefaultMessageLength)
        : source_file_name_(file)
        , times_overwritten_(0)   
        , source_file_line_(line) 
      {
          hash_ = combined_hash(file, line);
          message_ = new char[message_length];
#ifdef  DEBUG          
          if (1) printf("# WarningRecord:constructor allocates a new warning message string with max. %d chars"
              " at %p\n# ... for warnings launched at %s:%d --> hash = %16llx\n", 
          message_length, (void*)message_, file, line, hash_);
#endif
      } // constructor
      
      ~WarningRecord(void) {
#ifdef  DEBUG          
          if (1) printf("# WarningRecord:destructor: old warning message at %p for warnings launched at %s:%d reads:\n#\t%s\n", 
                            (void*)message_, get_sourcefile(), source_file_line_, message_);
#endif        
          // delete[] message_; // seems like this happens automagically
      } // destructor

      char* get_message(void) { ++times_overwritten_; return message_; }
      char* get_message_pointer(void) const { return message_; }
      char const* get_sourcefile(void) const { return source_file_name_.c_str(); }
      int get_sourceline(void) const { return source_file_line_; }
      size_t get_times(void) const { return times_overwritten_; }

  
  }; // class WarningRecord



  char* manage_warnings(char const *file, int const line, int const echo=0) {
    if (echo > 6) printf("\n# %s:%d  %s(file=%s, line=%d, echo=%d)\n", 
                      __FILE__, __LINE__, __func__, file, line, echo);

    static std::map<uint64_t,WarningRecord> map_;

    if (line < 1) { // line numbers created by the preprocessor start from 1

        assert('?' == file[0]); // make sure that we want special functionality
        if (0 == line) {
            // show_warnings has been called
            if (echo > 1) {
                auto const nw = map_.size();
                if ((echo < 3) || (nw < 1)) { // only give a summary of how many
                    printf("# %ld warnings have been recorded.\n", nw);
                } else {
                    printf("# %s: %ld recorded warnings:\n", __func__, nw);
                    size_t total_count = 0;
                    for (auto &hw : map_) {
                        auto const &w = hw.second;
                        auto const n_times = w.get_times();
                        printf("# \tin %s:%d (%ld times) \t%s\n", w.get_sourcefile(), 
                          w.get_sourceline(), n_times, w.get_message_pointer());
                        total_count += n_times;
                    } // w
                    if (nw > 0) printf("# %s: %ld warnings in total\n", __func__, total_count);
                } // summary
            } // echo
        } else {
            // clear_warnings has been called
            if (echo > 1) printf("# %s: clear all %ld entries from records\n", __func__, map_.size());
            map_.clear();
        }
        return nullptr;

    } else { // special functions

        // regular usage with file and line
        auto const hash = combined_hash(file, line);
        auto const search = map_.find(hash);
        if (map_.end() != search) {
            if (echo > 1) printf("# %s: found entry for hash %16lx\n", __func__, hash);
            return search->second.get_message();
        } else {
            if (echo > 1) printf("# %s: insert new entry for hash %16lx\n", __func__, hash);
            auto const iit = map_.insert({hash, WarningRecord(file, line)});
            // if (success) {
            return iit.first->second.get_message();
            // } else {
                // if (echo > 1) printf("# %s: failed, return new 256 chars\n", __func__);
                // return new char[256];
            // }
        } // found

    } // special functions

  } // manage_warnings

  char* new_warning(char const *file, int const line) {
      return manage_warnings(file, line);
  } // new_warning

  status_t show_warnings(int const echo) {
      return (nullptr != manage_warnings("?",  0, echo));
  } // show_warnings

  status_t clear_warnings(int const echo) {
      return (nullptr != manage_warnings("?", -1, echo));
  } // clear_warnings


#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_create_and_destroy(int const echo=9) {
    if (echo > 1) printf("\n# %s:%d  %s\n\n", __FILE__, __LINE__, __func__);
    WarningRecord wr(__FILE__,__LINE__);
    auto const msg = wr.get_message();
    sprintf(msg, "This is a non-recorded warning! Text created in %s:%d", __FILE__, __LINE__);
    return 0;
  } // test_create_and_destroy

  status_t test_preprocessor_macro(int const echo=9) {
    if (echo > 1) printf("\n# %s:%d  %s\n", __FILE__, __LINE__, __func__);
    warn("This is a test warning from %s:%d", __FILE__, __LINE__);
    return 0;
  } // test_preprocessor_macro

  status_t test_overwriting(int const echo=9) {
    if (echo > 1) printf("\n# %s:%d  %s\n", __FILE__, __LINE__, __func__);
    for(int i = 0; i < 9; ++i) {
        warn("This is a test warning from inside a loop, iteration #%d", i);
    } // i
    return 0;
  } // test_overwriting
  
  status_t all_tests() {
    auto status = 0;
    status += test_create_and_destroy();
    status += test_preprocessor_macro();
    status += test_overwriting();
    status += show_warnings(3); // display those warnings that have been launched for test purposes
    status += clear_warnings(2); // clear test warnings from record
    return status;
  } // all_tests
#endif // NO_UNIT_TESTS  

} // namespace recorded_warnings
