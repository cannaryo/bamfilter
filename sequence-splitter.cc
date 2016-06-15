// ---------------------------------------------------------------------------
// sequence-splitter.cc
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <sstream>
#include <iostream>
#include <vector>
#include "sequence-splitter.h"

namespace Constants = ::BamTools::Constants;
using BamTools::BamAlignment;
using BamTools::BamReader;
using BamTools::CigarOp;
using std::string;
using std::vector;

namespace fmu_tools {

// Open fastq file into fstream
bool SequenceSplitter::Open(string file) {
  output_.open( file.c_str(), std::ios::out );
  file_loaded_ = output_;
  return file_loaded_;
}

// Close fstream
void SequenceSplitter::Close() {
  if(output_) {
    output_.close();
    file_loaded_ = false;
  }
  return;
}

// Write data on file
bool SequenceSplitter::Write() {
  if(!file_loaded_ || !data_renewed_) {
    return false;
  }
  output_ << left_name_ << std::endl;
  output_ << left_sequence_ << std::endl;
  output_ << "+" << std::endl; 
  output_ << left_quality_ << std::endl;
  output_ << right_name_ << std::endl;
  output_ << right_sequence_ << std::endl;
  output_ << "+" << std::endl; 
  output_ << right_quality_ << std::endl;
  data_renewed_ = false;
  return true;
}

// Soft-clip split
bool SequenceSplitter::SplitBySoftClip(const BamTools::BamAlignment &d) {
  int len = d.Length;
  vector<CigarOp> cigar = d.CigarData;
  if(cigar.size() < 2) {
    return false;
  }
  int first = 0, last = cigar.size() - 1;
  int sz_l = 0, sz_r = 0;
  if(cigar[first].Type == 'S') {
    sz_l = cigar[first].Length;
  }
  if(cigar[last].Type == 'S') {
    sz_r = cigar[last].Length;
  }
  if(sz_l == 0 && sz_r == 0) {
    return false;
  }
  if(sz_l > sz_r) {
    if(sz_l < min_clip_ || (len - sz_l) < min_hold_) {
      return false;
    }
    left_sequence_ = d.QueryBases.substr(0, sz_l);
    left_quality_ = d.Qualities.substr(0, sz_l);
    left_name_ = d.Name + ":LS";
    right_sequence_ = d.QueryBases.substr(sz_l, len - sz_l);
    right_quality_ = d.Qualities.substr(sz_l, len - sz_l);
    right_name_ = d.Name + ":RM";
  } else {
    if(sz_r < min_clip_ || (len - sz_r) < min_hold_) {
      return false;
    }
    left_sequence_ = d.QueryBases.substr(0, len - sz_r);
    left_quality_ = d.Qualities.substr(0, len - sz_r);
    left_name_ = d.Name + ":LM";
    right_sequence_ = d.QueryBases.substr(len - sz_r, sz_r);
    right_quality_ = d.Qualities.substr(len - sz_r, sz_r);
    right_name_ = d.Name + ":RS";
  }
  if(d.Name[0] != '@') {
    right_name_ = "@" + right_name_;
    left_name_ = "@" + left_name_;
  }
  data_renewed_ = true;
  return true;
}

bool SequenceSplitter::SplitByFixedLength(const BamTools::BamAlignment &d, int size) {
  int len = d.Length;
  if(len < min_length_ || len < size*2) {
    return false;
  }
  left_sequence_ = d.QueryBases.substr(0, size);
  left_quality_ = d.Qualities.substr(0, size);
  left_name_ = d.Name + ":LF";
  right_sequence_ = d.QueryBases.substr(len - size, size);
  right_quality_ = d.Qualities.substr(len -size, size);
  right_name_ = d.Name + ":RM";

  if(d.Name[0] != '@') {
    right_name_ = "@" + right_name_;
    left_name_ = "@" + left_name_;
  }
  data_renewed_ = true;
  return true;
}

} // fmu_tools
