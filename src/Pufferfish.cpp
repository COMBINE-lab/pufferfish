//
// Pufferfish - An efficient compacted dBG index
// 
// Copyright (C) 2017 Rob Patro, Fatemeh Almodaresi, Hirak Sarkar
//
// This file is part of Pufferfish.
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

#include <iostream>
#include <fstream>
#include <vector>

//#include <cereal/archives/json.hpp>

#include "PufferfishConfig.hpp"
//#include "IndexHeader.hpp"

int pufferfishIndex(int argc, char* argv[]);
int pufferfishTest(int argc, char* argv[]);
//int rapMapMap(int argc, char* argv[]);
//int rapMapSAMap(int argc, char* argv[]);

void printUsage() {
    std::string versionString = pufferfish::version;
    std::cerr << "Pufferfish v"
              << versionString << '\n';
    std::cerr << "=====================================\n";
    auto usage =
        R"(
There are currently 2 Pufferfish subcommands
    index --- builds a Pufferfish index
    test --- tests k-mer lookup in the index

Run a corresponding command "pufferrish <cmd> -h" for
more information on each of the possible Pufferfish
commands.)";
    std::cerr << usage << '\n';
}


int main(int argc, char* argv[]) {

    std::vector<char*> args;
    args.push_back(argv[0]);

    if (argc < 2) {
        printUsage();
        std::exit(0);
    }

    for (int i = 2; i < argc; ++i) {
        args.push_back(argv[i]);
    }

    if (std::string(argv[1]) == "-h" or
        std::string(argv[1]) == "--help") {
        printUsage();
        std::exit(0);
    }

    if (std::string(argv[1]) == "index") {
        return pufferfishIndex(argc - 1, args.data());
    } else if (std::string(argv[1]) == "test") {
        return pufferfishTest(argc - 1, args.data());
    } else {
        std::cerr << "the command " << argv[1]
                  << " is not yet implemented\n";
      return 1;
    }
    return 0;
}
