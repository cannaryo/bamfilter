// ---------------------------------------------------------------------------
// sequence-splitter.h
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------
// Split tool for sequence data
// ---------------------------------------------------------------------------

#ifndef FMUTOOLS_SEQUENCESPLITTER_H_
#define FMUTOOLS_SEQUENCESPLITTER_H_

#include <string>
#include <fstream>
#include <api/BamConstants.h>
#include <api/BamReader.h>
#include <api/BamAlignment.h>

namespace fmu_tools {

// This class supplies sequence-split operations for BAM/SAM data.
// The split sequences can be converted into FASTQ format
class SequenceSplitter {
 public:
  SequenceSplitter() {
    min_clip_ = 0;
    min_hold_ = 0;
    min_length_ = 0;
    file_loaded_ = false;
    data_renewed_ = false;
  }
  ~SequenceSplitter() {}

 public:
  // Opens Fastq file.
  // When it failed to open the file, returns false;
  bool Open(std::string file);
  // Closes Fastq file if it is opend.
  void Close();
  // Writes present data on file.
  // When the data is not renewed or the file is not opened, it returns false.  
  bool Write();
  // Makes soft-clip split sequence from BAM alignment data.
  // If it successfully split the sequence, data member are renewed and returns true.  
  // If data do not satisfy the condition, just returns false.
  bool SplitBySoftClip(const BamTools::BamAlignment &d);
  // Makes Fixed-length split sequence from BAM alignment data.
  // If it successfully split the sequence, data member are renewed and returns true.  
  // If data do not satisfy the condition, just returns false.
  bool SplitByFixedLength(const BamTools::BamAlignment &d);

  // Accessors
  void set_min_clip(int n) { min_clip_ = n; }
  void set_min_hold(int n) { min_hold_ = n; }
  void set_min_length(int n) { min_length_ = n; }
  
 private:
  int min_clip_;
  int min_hold_;
  int min_length_;
  bool file_loaded_;
  bool data_renewed_; 
  std::ofstream output_;
  std::string left_name_;
  std::string left_sequence_;
  std::string left_quality_;
  std::string right_name_;
  std::string right_sequence_;
  std::string right_quality_;
};

} // fmu_tools

#endif  // FMUTOOLS_SEQUENCESPLITTER_H_
