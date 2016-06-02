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
#include <api/BamConstants.h>
#include <api/BamReader.h>
#include <api/BamAlignment.h>

namespace fmu_tools {

// This class supplies sequence-split operations for BAM/SAM data.
// The split sequences can be converted into FASTQ format
class SequenceSplitter {
 public:
  SequenceSplitter() {}
  ~SequenceSplitter() {}

 public:
  // Set reference id and header values.
  // Read access (d) to BAM file must be specified.
  //  void SetFileInformation(const BamTools::BamReader& d);
  // Converts BAM alignment (d) into SAM format and returns it as string. 
  //  std::string ConvertToSam(const BamTools::BamAlignment& d);

 private:
  //  BamTools::RefVector reference_names_;
  //  BamTools::SamHeader sam_header_;

};

} // fmu_tools

#endif  // FMUTOOLS_SEQUENCESPLITTER_H_
