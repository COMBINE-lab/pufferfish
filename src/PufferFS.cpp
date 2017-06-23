//
// RapMap - Rapid and accurate mapping of short reads to transcriptomes using
// quasi-mapping.
// Copyright (C) 2015, 2016 Rob Patro, Avi Srivastava, Hirak Sarkar
//
// This file is part of RapMap.
//
// RapMap is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// RapMap is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with RapMap.  If not, see <http://www.gnu.org/licenses/>.
//

#include "PufferFS.hpp"
#include <sys/stat.h>

namespace puffer {
namespace fs {

// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool FileExists(const char* path) {
  struct stat fileStat;
  if (stat(path, &fileStat)) {
    return false;
  }
  if (!S_ISREG(fileStat.st_mode)) {
    return false;
  }
  return true;
}

// Taken from
// http://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool DirExists(const char* path) {
  struct stat fileStat;
  if (stat(path, &fileStat)) {
    return false;
  }
  if (!S_ISDIR(fileStat.st_mode)) {
    return false;
  }
  return true;
}

void MakeDir(const char* path) { mkdir(path, ACCESSPERMS); }
}
}
