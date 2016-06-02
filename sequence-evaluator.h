// ---------------------------------------------------------------------------
// sequence-filter.h
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------
// Filter tool for sequence data
// ---------------------------------------------------------------------------

#ifndef FMUTOOLS_SEQUENCEEVALUATOR_H_
#define FMUTOOLS_SEQUENCEEVALUATOR_H_

#include <string>
#include <api/BamConstants.h>
#include <api/BamReader.h>
#include <api/BamAlignment.h>

namespace fmu_tools {

// This class contains various operations for sequece data.
// The supporting input is BAM alignment format
// See indivisual comment before each method.   
class SequenceEvaluator {
 public:
  SequenceEvaluator() {}
  ~SequenceEvaluator() {}

 public:
  struct FilterParameters {
    int min_del; 
    int min_ins; 
    int min_indel; 
    int min_softclip; 
    int min_match;
    bool keep_unmapped; 
    bool reverse_sign;
  };

  // Set reference id and header values.
  // Read access (d) to BAM file must be specified.
  void SetFileInformation(const BamTools::BamReader &d);

  // Converts BAM alignment (d) into SAM format and returns it as string. 
  std::string ConvertToSam(const BamTools::BamAlignment &d) const;

  // Returns true if the data satisfy the conditions and should be kept.
  // It works as a filter of BAM Alignment data by CIGAR count.
  // Using BamReader::GetNextAlignmentCore for input data 
  //   instead of BamReader::GetNextAlignment is effective for speed.
  bool FilterByCigar(const BamTools::BamAlignment &d, const FilterParameters &prms) const;

 private:
  BamTools::RefVector reference_names_;
  BamTools::SamHeader sam_header_;

};

} // fmu_tools

#endif  // FMUTOOLS_SEQUENCEEVALUATOR_H_
