// ---------------------------------------------------------------------------
// bamfilter.h
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------
// Fast filter for BAM format file
// ---------------------------------------------------------------------------

#ifndef FMUTOOLS_BAMFILTER_H_
#define FMUTOOLS_BAMFILTER_H_

#include <string>

namespace fmu_tools {

// Main framework of bamfilter 
class BamFilterMain {

 public:
  BamFilterMain() {}
  ~BamFilterMain() {}
  
 public:
  // Main entry point
  int Run(int argc, char* argv[]);

};

} // fmu_tools

#endif  // FMUTOOLS_BAMFILTER_H_
