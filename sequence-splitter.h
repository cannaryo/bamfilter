// ---------------------------------------------------------------------------
// sequence-splitterh
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------
// Utilities for sequence data
// ---------------------------------------------------------------------------

#ifndef FMUTOOLS_SEQUENCESPLITTER_H_
#define FMUTOOLS_SEQUENCESPLITTER_H_

#include <string>

//extern class BamTools::BamAlignment;

namespace fmu_tools {

// This class contains various operations for sequece data.
// See indivisual comment before each method.   
class SequenceSplitter {
 public:
  SequenceSplitter() {}
  ~SequenceSplitter() {}

 public:
  // Returns string as general SAM format. 
  std::string ConvertToSam(const BamTools::BamAlignment& d);

};
  
} // fmu_tools

#endif  // FMUTOOLS_SEQUENCESPLITTER_H_
