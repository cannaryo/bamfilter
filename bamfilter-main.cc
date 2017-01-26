// ---------------------------------------------------------------------------
// bamfilter.cc
// Copyright (c) 2016 - Ryo Kanno
// ml_kanno@csc.jp
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
// ---------------------------------------------------------------------------

#include <stdio.h>
#include <getopt.h>
#include <sstream>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <cmdline.h>
#include <api/BamReader.h>
#include "sequence-splitter.h"
#include "sequence-evaluator.h"
#include "bamfilter-main.h"

namespace Constants = ::BamTools::Constants;
using BamTools::BamReader;
using BamTools::BamAlignment;
using BamTools::CigarOp;
using std::string;
using std::vector;

namespace fmu_tools {

// ---------------------------------------------
// OptionHandler implementation

class BamFilterMain::OptionHandler {
public:
  OptionHandler() {};
  ~OptionHandler() {};

public:
  bool ParseOption(int argc, char *argv[]);
  bool IsOption(const string &name) {
    return cmd_parser_.exist(name);
  }
  int GetInt(const string &name) {
    return cmd_parser_.get<int>(name);
  }
  string GetString(const string &name) {
    return cmd_parser_.get<string>(name);
  }
  string input_file() {
    return input_file_;
  }

private:
  void SetOptionArg();

  cmdline::parser cmd_parser_;
  string input_file_;
};

// Parse command line options and initialize parameters
bool BamFilterMain::OptionHandler::ParseOption(int argc, char *argv[]) {
  this->SetOptionArg();
  if( !cmd_parser_.parse(argc, argv) ) {
    std::cerr << cmd_parser_.error_full() << std::endl;
    std::cerr << "See '--help' for valid options." << std::endl;
    return false;
  }
  vector<string> rest = cmd_parser_.rest();
  if( this->IsOption("help") || rest.size() < 1 ) {
    std::cerr << cmd_parser_.usage() << std::endl;
    return false;
  }
  input_file_ = rest[0];
  return true;
}

// Set all options for the application
void BamFilterMain::OptionHandler::SetOptionArg()
{
  cmd_parser_.footer("<in.bam>: Filter and split read data in BAM file");
  cmd_parser_.add("help", 0, "print help message");
  cmd_parser_.add<int>("min-deletion", 'd', "keep deletion >= N", false, -1);
  cmd_parser_.add<int>("min-insertion", 'i', "keep insertion >= N", false, -1);
  cmd_parser_.add<int>("min-indel", 'b', "keep indel >= N", false, -1);
  cmd_parser_.add<int>("min-softclip", 's', "keep softclip >= N", false, -1);
  cmd_parser_.add<int>("min-match", 'm', "keep match >= N", false, -1);
  cmd_parser_.add("reverse", 'r', "reverse inequality sign for -[dibsm]");
  cmd_parser_.add("unmapped",'u', "keep unmapped reads");
  cmd_parser_.add("primary", 'p', "only primary alignment");
  cmd_parser_.add<string>("softclip", 0, "apply soft-clip split (output: FASTQ file)", false);
  cmd_parser_.add<int>("soft-min-clip", 'c', "minimum length of clipped side (default: 25)", false, 25);
  cmd_parser_.add<int>("soft-min-hold", 'h', "minimum length of remaining side (default: 25)", false, 25);
  cmd_parser_.add<string>("fixed", 0, "apply fixed-lenght split (output: FASTQ file)", false);
  cmd_parser_.add<int>("fix-min-length", 'l', "minimum read length for fixed-length split (default: 100)", false, 100);
  cmd_parser_.add<int>("fix-split-size", 'f', "split size (default: 35)", false, 35);
  cmd_parser_.add("keep-only-final", 'k', "do not output intermediate reads (it saves time)");
  cmd_parser_.add<string>("output", 'o', "output filtered data in SAM format (output: SAM file)", false);
}

// ---------------------------------------------
// BamFilterMain implementation

int BamFilterMain::Run(int argc, char* argv[]) {
  if( !opt().ParseOption(argc, argv) ) {
    return 0;
  }
  BamReader file;
  if( !file.Open(opt().input_file()) ) {
    std::cerr << "Fail to open file: " << opt().input_file() << std::endl;
    return 0;
  }
  std::cerr << "Input: " <<  file.GetFilename() << std::endl;
  
  BamAlignment record;
  std::ofstream output;
  string head = file.GetHeaderText();
  SequenceEvaluator seq_filter;
  struct SequenceEvaluator::FilterParameters prms;
  prms.min_del = opt().GetInt("min-deletion");
  prms.min_ins = opt().GetInt("min-insertion");
  prms.min_indel = opt().GetInt("min-indel");
  prms.min_softclip = opt().GetInt("min-softclip");
  prms.min_match = opt().GetInt("min-match");
  prms.reverse_sign = opt().IsOption("reverse");
  prms.keep_unmapped = opt().IsOption("unmapped");
  SequenceSplitter soft_split;
  bool do_softclip = opt().IsOption("softclip");
  if( do_softclip ) {
    if( !soft_split.Open(opt().GetString("softclip")) ) {
      return 1;
    }
    soft_split.set_min_clip( opt().GetInt("soft-min-clip") );
    soft_split.set_min_hold( opt().GetInt("soft-min-hold") );
  }
  SequenceSplitter fix_split;
  int fix_size = opt().GetInt("fix-split-size");
  bool do_fixed = opt().IsOption("fixed");
  if( do_fixed ) {
    if ( !fix_split.Open(opt().GetString("fixed")) ) {
      return 1;
    }
    fix_split.set_min_length( opt().GetInt("fix-min-length") );
  }
  bool do_out = opt().IsOption("output");
  if( do_out ) {
    output.open( opt().GetString("output").c_str(), std::ios::out );
    output << head;
  }
  bool keep_only_final = opt().IsOption("keep-only-final");
  bool primary_only = opt().IsOption("primary");
  seq_filter.SetFileInformation(file);

  int cc_all=0, cc_filt=0, cc_sc=0, cc_um=0, cc_out=0;
  while( file.GetNextAlignmentCore(record) ) {
    ++cc_all;
    if( primary_only && !record.IsPrimaryAlignment() ) {
      continue;
    }
    if( seq_filter.FilterByCigar(record, prms) ) {
      ++cc_filt;
      record.BuildCharData();

      if( do_softclip && soft_split.SplitBySoftClip(record) ) {
        ++cc_sc;
        soft_split.Write();
        if( keep_only_final && do_fixed ) {
          continue;
        }
      } else if( do_fixed && fix_split.SplitByFixedLength(record, fix_size) ) {
        ++cc_um;
        fix_split.Write();
      } else if( keep_only_final && (do_fixed || do_softclip) ) {
        continue;
      }
      if( do_out ) {
        // write alignment on SAM file
        string str = seq_filter.ConvertToSam(record);
        ++cc_out;
        output << str;
      }
    }
  } // while loop

  std::cerr << "Keep records: " << cc_filt << " / " << cc_all << std::endl;

  if( opt().IsOption("softclip") ) {
    std::cerr << cc_sc << " records were written in " << opt().GetString("softclip") << std::endl;
    soft_split.Close();
  }
  if( opt().IsOption("fixed") ) {
    std::cerr << cc_um << " records were written in " << opt().GetString("fixed") << std::endl;    
    fix_split.Close();
  }
  if( opt().IsOption("output") ) {
    std::cerr << cc_out << " records were written in " << opt().GetString("output") << std::endl;    
  }
  file.Close();
  
  return 0;
}

// Constructor
BamFilterMain::BamFilterMain() {
  opt_handler_ = new OptionHandler();
}
 
// Destructor
BamFilterMain::~BamFilterMain() {
  delete opt_handler_;
  opt_handler_ = NULL;
}

} // fmu_tools
