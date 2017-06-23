#include "Util.hpp"

namespace util {

char complement(char& c) {
  switch (c) {
  case 'A':
    c = 'T';
    return c;
  case 'T':
    c = 'A';
    return c;
  case 'C':
    c = 'G';
    return c;
  case 'G':
    c = 'C';
    return c;
  }
  return 'A';
}

std::string revcomp(std::string s) {
  int n = s.size();
  int halfLength = s.size() / 2;
  for (int i = 0; i < halfLength; i++) {
    char temp = complement(s[i]);
    s[i] = complement(s[n - 1 - i]);
    s[n - 1 - i] = temp;
  }
  if (s.size() % 2 != 0) {
    s[halfLength] = complement(s[halfLength]);
  }
  return s;
}

bool isRevcomp(std::string s) {
  int n = s.size();
  int halfLength = n / 2;
  for (int i = 0; i < halfLength; i++) {
    char temp = complement(s[n - 1 - i]);
    if (temp != s[i])
      return false;
  }
  return true;
}
std::vector<std::pair<std::string, bool>> explode(const stx::string_view str,
                                                  const char& ch) {
  std::string next;
  std::vector<std::pair<std::string, bool>> result;
  // For each character in the string
  for (auto it = str.begin(); it != str.end(); it++) {
    // If we've hit the terminal character
    if (*it == '+' or *it == '-') {
      bool orientation = true;
      // If we have some characters accumulated
      // Add them to the result vector
      if (!next.empty()) {
        if (*it == '-') {
          orientation = false;
        }
        result.emplace_back(next, orientation);
        next.clear();
      }
    } else if (*it != ch) {
      // Accumulate the next character into the sequence
      next += *it;
    }
  }
  if (!next.empty())
    result.emplace_back(next, true); // this case shouldn't even happen
  return result;
}

bool is_number(const std::string& s) {
  return !s.empty() && std::find_if(s.begin(), s.end(), [](char c) {
                         return !std::isdigit(c);
                       }) == s.end();
}

// Avoiding un-necessary stream creation + replacing strings with string view
// is a bit > than a 2x win!
// implementation from : https://marcoarena.wordpress.com/tag/string_view/
std::vector<stx::string_view> split(stx::string_view str, char delims) {
  std::vector<stx::string_view> ret;

  stx::string_view::size_type start = 0;
  auto pos = str.find_first_of(delims, start);
  while (pos != stx::string_view::npos) {
    if (pos != start) {
      ret.push_back(str.substr(start, pos - start));
    }
    start = pos + 1;
    pos = str.find_first_of(delims, start);
  }
  if (start < str.length()) {
    ret.push_back(str.substr(start, str.length() - start));
  }
  return ret;
}
}
