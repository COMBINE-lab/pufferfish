#ifndef _UTIL__H
#define _UTIL__H

#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include "string_view.hpp"


namespace util {

  class IndexOptions {
  public:
    uint32_t k{31};
    std::string gfa_file;
    std::string cfile;
    std::string rfile;
    std::string outdir;
  };

  class TestOptions {
  public:
  };

  char complement(char& c);
  std::string revcomp(std::string s);
  std::vector<std::pair<std::string, bool> > explode(const stx::string_view str, const char& ch);
	bool is_number(const std::string& s);
  // Avoiding un-necessary stream creation + replacing strings with string view
  // is a bit > than a 2x win!
  // implementation from : https://marcoarena.wordpress.com/tag/string_view/
  std::vector<stx::string_view> split(stx::string_view str, char delims);
}

#endif // _UTIL__H
