//
// Created by Andrey Bzikadze on 2/20/21.
//

#pragma once

#include <fstream>
#include <string>

inline std::ostream &operator<<(std::ostream &os, unsigned __int128 val) {
  std::vector<size_t> res;
  while (val != 0) {
    res.push_back(val % 10);
    val /= 10;
  }
  for (auto it = res.rbegin(); it != res.rend(); ++it) { os << *it; }
  return os;
}

inline std::istream &operator>>(std::istream &is, unsigned __int128 &val) {
  val = 0;
  std::string s;
  is >> s;
  for (char &it : s) {
    val *= 10;
    val += it - '0';
  }
  return is;
}
