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
  bool Exist(const string &name) {
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
  if(!cmd_parser_.parse(argc, argv)) {
    std::cerr << cmd_parser_.error_full() << std::endl;
    std::cerr << "See '--help' for valid options." << std::endl;
    return false;
  }
  vector<string> rest = cmd_parser_.rest();
  if(this->Exist("help") || rest.size() < 1) {
    std::cerr << cmd_parser_.usage() << std::endl;
    return false;
  }
  input_file_ = rest[0];
  return true;
}

// Set all options for the application
void BamFilterMain::OptionHandler::SetOptionArg()
{
  cmd_parser_.footer("<in.bam>: Filter and split BAM file");
  cmd_parser_.add("help", 0, "print help message");
  cmd_parser_.add<int>("min-deletion", 'd', "keep deletion >= N", false, -1);
  cmd_parser_.add<int>("min-insertion", 'i', "keep insertion >= N", false, -1);
  cmd_parser_.add<int>("min-indel", 'b', "keep indel >= N", false, -1);
  cmd_parser_.add<int>("min-softclip", 's', "keep softclip >= N", false, -1);
  cmd_parser_.add<int>("min-match", 'm', "keep match >= N", false, -1);
  cmd_parser_.add("reverse", 'r', "reverse inequality sign for -[dibsm]");
  cmd_parser_.add("unmapped",'u', "keep unmapped reads");
  cmd_parser_.add<string>("softclip", 0, "apply soft-clip split (output: FASTQ file)", false);
  cmd_parser_.add<int>("soft-min-clip", 'c', "minimum length of clipped side (default: 25)", false, 25);
  cmd_parser_.add<int>("soft-min-hold", 'h', "minimum length of remaining side (default: 25)", false, 25);
  cmd_parser_.add("keep-used-read", 'k', "keep used reads for soft-clip split (generally not necessary)");
  cmd_parser_.add<string>("fixed", 0, "apply fixed-lenght split (output: FASTQ file)", false);
  cmd_parser_.add<int>("fix-min-length", 'l', "minimum read length for fixed-length split (default: 100)", false, 100);
  cmd_parser_.add<int>("fix-split-size", 'p', "split size (default: 35)", false, 35);
  cmd_parser_.add<string>("output", 'o', "output filtered data (output: SAM file)", false);
}

// ---------------------------------------------
// BamFilterMain implementation

int BamFilterMain::Run(int argc, char* argv[]) {
  if(!opt().ParseOption(argc, argv)) {
    return 0;
  }
  BamReader file;
  if(!file.Open( opt().input_file() )) {
    std::cerr << "Fail to open file: " << opt().input_file() << std::endl;
    return 0;
  }
  std::cerr << "Input: " <<  file.GetFilename() << std::endl;
  
  BamAlignment record;
  std::ofstream output;
  bool do_out = false;
  string head = file.GetHeaderText();
  SequenceEvaluator seq_filter;
  struct SequenceEvaluator::FilterParameters prms;
  prms.min_del = opt().GetInt("min-deletion");
  prms.min_ins = opt().GetInt("min-insertion");
  prms.min_indel = opt().GetInt("min-indel");
  prms.min_softclip = opt().GetInt("min-softclip");
  prms.min_match = opt().GetInt("min-match");
  prms.reverse_sign = opt().Exist("reverse");
  prms.keep_unmapped = opt().Exist("unmapped");

  if( opt().Exist("output") ) {
    output.open( opt().GetString("output").c_str(), std::ios::out );
    do_out = true;
    output << head;
  }
  seq_filter.SetFileInformation(file);
  while(file.GetNextAlignmentCore(record)) {
    if(seq_filter.FilterByCigar(record, prms)) {
      record.BuildCharData();
      string str = seq_filter.ConvertToSam(record);
      if(do_out) {
        output << str;
      } else {
        std::cout << str;
      }
    }
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
