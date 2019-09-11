#pragma once

#include "chemical_symbol.h" // element_symbols

typedef int status_t;

namespace chemical_symbol {
  
// A    --> l:13, r:18, s:33, g:47, u:79, t:85, c:89, m:95, 
// B    --> e:4,  :5, r:35, a:56, i:83, k:97, h:107, 
// C    -->  :6, l:17, a:20, r:24, o:27, u:29, d:48, s:55, e:58, m:96, f:98, n:112, 
// D    --> y:66, b:105, s:110, 
// E    --> u:63, r:68, s:99, 
// F    -->  :9, e:26, r:87, m:100, 
// G    --> a:31, e:32, d:64, 
// H    -->  :1, e:2, o:67, f:72, g:80, s:108, 
// I    --> n:49,  :53, r:77, 
// J    --> 
// K    -->  :19, r:36, 
// L    --> i:3, a:57, u:71, r:103, 
// M    --> g:12, n:25, o:42, d:101, t:109, 
// N    -->  :7, e:10, a:11, i:28, b:41, d:60, p:93, o:102, 
// O    -->  :8, s:76, 
// P    -->  :15, d:46, r:59, m:61, t:78, b:82, o:84, a:91, u:94, 
// Q    --> 
// R    --> b:37, u:44, h:45, e:75, n:86, a:88, f:104, g:111, 
// S    --> i:14,  :16, c:21, e:34, r:38, n:50, b:51, m:62, g:106, 
// T    --> i:22, c:43, e:52, b:65, m:69, a:73, l:81, h:90, 
// U    -->  :92, 
// V    -->  :23, 
// W    -->  :74, 
// X    --> e:54, 
// Y    -->  :39, b:70, 
// Z    --> n:30, r:40, 
// _    --> _:0, 
// e    -->  :127, 
// u    --> t:113, q:114, p:115, h:116, s:117, o:118, n:119, d:120, 1:121, 2:122, 3:123, 4:124, 5:125, 6:126, 
  
//   inline int decode(char const S, char const y) {
//       switch (S) { // first character, usually upper case
//                case 'A':
//         break; case 'B':
//         break; case 'C':
//         break; case 'D':
//         break; case 'E':
//         break; case 'F':
//         break; case 'G':
//         break; case 'H':
//         break; case 'I':
//         break; case 'J':
//         break; case 'K':
//         break; case 'L':
//         break; case 'M':
//         break; case 'N':
//         break; case 'O':
//         break; case 'P':
//         break; case 'Q':
//         break; case 'R':
//         break; case 'S':
//         break; case 'T':
//         break; case 'U':
//         break; case 'V':
//         break; case 'W':
//         break; case 'X':
//         break; case 'Y':
//         break; case 'Z':
//         break; case '_':
//         break; case 'e':
//         break; case 'u':
//         break;  default:
//       } // S
//       return 0;
//   } // decode
  
#ifdef  NO_UNIT_TESTS
  status_t all_tests() { printf("\nError: %s was compiled with -D NO_UNIT_TESTS\n\n", __FILE__); return -1; }
#else // NO_UNIT_TESTS

  status_t test_decoding(int const echo=5) {
    
      int8_t Z[128][16];
      uint8_t nZ[128];
      
      for(int c = 0; c < 128; ++c) {
          nZ[c] = 0;
          for(int n = 0; n < 16; ++n) Z[c][n] = -1;
      } // c

      for(int m = 0; m < 128; ++m) {
          int c = element_symbols[2*m + 0];
          Z[c][nZ[c]] = m;
          ++nZ[c];
      } // m
      
      for(int c = 65; c < 118; ++c) {
          printf("# \'%c\' --> ", c);
          for(int n = 0; n < nZ[c]; ++n) {
              int const m = Z[c][n];
              printf("\'%c\':%d, ", element_symbols[2*m + 1], m);
          } // n
          printf("\n");
      } // c
      
      return 0;
  } // test_decoding

  status_t test_decoding_char2int(int const echo=5) {
      
      for(int Z = 0; Z < 128; ++Z) {
          printf("        case %d: return %d;\n", (element_symbols[2*Z + 0] - 48) * 128 + (element_symbols[2*Z + 1] - 48), Z);
      } // Z
      
      return 0;
  } // test_decoding
  
  status_t all_tests() {
    auto status = 0;
    status += test_decoding_char2int();
    return status;
  } // all_tests

#endif // NO_UNIT_TESTS

} // namespace chemical_symbol
