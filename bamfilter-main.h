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
#include <memory>

namespace fmu_tools {

// Main framework of bamfilter 
class BamFilterMain {

 public:
  BamFilterMain(); 
  ~BamFilterMain();
  
 public:
  // Main entry point
  int Run(int argc, char *argv[]);

 private:
  // Access to command line options
  // Use only for clearly defined operations
  class OptionHandler;
  OptionHandler *opt_handler_;
  OptionHandler &opt() { return *opt_handler_; }

};

} // fmu_tools

#endif  // FMUTOOLS_BAMFILTER_H_
