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
#include <sstream>
#include <iostream>
#include <string>
#include <vector>
#include <api/BamReader.h>
#include "sequence-splitter.h"
#include "bamfilter-main.h"

namespace Constants = ::BamTools::Constants;
using BamTools::BamReader;
using BamTools::BamAlignment;
using BamTools::CigarOp;
using std::string;
using std::vector;

namespace fmu_tools {

int BamFilterMain::Run(int argc, char* argv[]) {
  BamReader* file = new BamTools::BamReader();
  BamAlignment ali;
  file->Open("testdata/s1_test_1_1_20.bam");
  printf("File %s is opend\n", file->GetFilename().c_str());  
  string head = file->GetHeaderText();
  printf("%s\n",head.c_str());
  int i = 0;
  SequenceSplitter seq;
  while(file->GetNextAlignment(ali)) {
    string str = seq.ConvertToSam(ali);
    printf("%s\n", str.c_str());
    if(++i > 10)
      break;
  }
  file->Close();
  delete file;
  return 0;
}

} // fmu_tools
