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
#include <string>
#include <vector>
#include "sequence-evaluator.h"

namespace Constants = ::BamTools::Constants;
using BamTools::BamAlignment;
using BamTools::BamReader;
using BamTools::CigarOp;
using std::string;
using std::vector;

namespace {

void CountCigar(const vector<CigarOp> &cigars, vector<int> *cigar_count) {
  // Order of CIGAR count
  // 0:M, 1:I, 2:D, 3:N, 4:S, 5:H, 6:P, 7:X, 8:=
  for ( vector<CigarOp>::const_iterator iter=cigars.begin(); iter != cigars.end(); ++iter ) {
    switch(iter->Type) {
    case 'M':
      (*cigar_count)[0] += iter->Length;
      break;
    case 'I':
      (*cigar_count)[1] += iter->Length;
      break;
    case 'D':
      (*cigar_count)[2] += iter->Length;
      break;
    case 'N':
      (*cigar_count)[3] += iter->Length;
      break;
    case 'S':
      (*cigar_count)[4] += iter->Length;
      break;
    case 'H':
      (*cigar_count)[5] += iter->Length;
      break;
    case 'P':
      (*cigar_count)[6] += iter->Length;
      break;
    case 'X':
      (*cigar_count)[7] += iter->Length;
      break;
    case '=':
      (*cigar_count)[8] += iter->Length;
      break;      
    }
  }
}

} // namespace

namespace fmu_tools {

// Filter by CIGAR count
bool SequenceEvaluator::FilterByCigar(const BamAlignment &d, const FilterParameters &prms) const {
  // Order of CIGAR count
  // 0:M, 1:I, 2:D, 3:N, 4:S, 5:H, 6:P, 7:X, 8:=
  vector<int> ncigars(9,0);
  CountCigar(d.CigarData, &ncigars);
  if(prms.keep_unmapped && !d.IsMapped()) {
    return true;
  } else if(prms.reverse_sign) {
    if(prms.min_match != -1 && prms.min_match >= ncigars[0]) {
      return true;
    } else if(prms.min_ins != -1 && prms.min_ins >= ncigars[1]) {
      return true;
    } else if(prms.min_del != -1 && prms.min_del >= ncigars[2]) {
      return true;
    } else if(prms.min_indel != -1 && prms.min_indel >= ncigars[1] + ncigars[2]) {
      return true;
    } else if(prms.min_softclip != -1 && prms.min_softclip >= ncigars[4]) {
      return true;
    }
  } else {
    if(prms.min_match != -1 && prms.min_match <= ncigars[0]) {
      return true;
    } else if(prms.min_ins != -1 && prms.min_ins <= ncigars[1]) {
      return true;
    } else if(prms.min_del != -1 && prms.min_del <= ncigars[2]) {
      return true;
    } else if(prms.min_indel != -1 && prms.min_indel <= ncigars[1] + ncigars[2]) {
      return true;
    } else if(prms.min_softclip != -1 && prms.min_softclip <= ncigars[4]) {
      return true;
    }
  }
  return false;  
}

// Attach reference name and header information to the object
void SequenceEvaluator::SetFileInformation(const BamReader &d) {
  reference_names_ = d.GetReferenceData();
  sam_header_ = d.GetHeader();
}

// print BamAlignment in SAM format
string SequenceEvaluator::ConvertToSam(const BamAlignment &d) const {
  // <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [ <TAG>:<VTYPE>:<VALUE> [...] ]
  // write name & alignment flag
  std::ostringstream out;
  out << d.Name << "\t" << d.AlignmentFlag << "\t";
  
  // write reference name
  if ( (d.RefID >= 0) && (d.RefID < (int)reference_names_.size()) ) 
    out << reference_names_[d.RefID].RefName << "\t";
  else 
    out << "*\t";
  
  // write position & map quality
  out << d.Position+1 << "\t" << d.MapQuality << "\t";
    
  // write CIGAR
  const vector<CigarOp>& cigars = d.CigarData;
  if ( cigars.empty() ) 
    out << "*\t";
  else {
    for (vector<CigarOp>::const_iterator iter=cigars.begin() ; iter != cigars.end(); ++iter ) {
      //      const CigarOp& op = (*iter);
      out << iter->Length << iter->Type;
    }
    out << "\t";
  }

  // write mate reference name, mate position, & insert size
  if ( d.IsPaired() && (d.MateRefID >= 0) && (d.MateRefID < (int)reference_names_.size()) ) {
    if ( d.MateRefID == d.RefID )
      out << "=\t";
    else
      out << reference_names_[d.MateRefID].RefName << "\t";
    out << d.MatePosition+1 << "\t" << d.InsertSize << "\t";
  } 
  else
    out << "*\t0\t0\t";
    
  // write sequence
  if ( d.QueryBases.empty() )
    out << "*\t";
  else
    out << d.QueryBases << "\t";
    
  // write qualities
  if ( d.Qualities.empty() || (d.Qualities.at(0) == (char)0xFF) )
    out << "*";
  else
    out << d.Qualities;
    
  // write tag data
  const char* tag_data = d.TagData.c_str();
  const size_t tag_length = d.TagData.length();
  size_t index = 0;

  while ( index < tag_length ) {
    // write tag name   
    string tag_name = d.TagData.substr(index, 2);
    out << "\t" << tag_name << ":";
    index += 2;        
    // get data type    
    char type = d.TagData.at(index);
    ++index;
    switch ( type ) {
    case (Constants::BAM_TAG_TYPE_ASCII) :
      out << "A:" << tag_data[index];
      ++index;
      break;      
    case (Constants::BAM_TAG_TYPE_INT8) :
      // force value into integer-type (instead of char value)
      out << "i:" << static_cast<int16_t>(static_cast<int8_t>(tag_data[index]));
      ++index;
      break;      
    case (Constants::BAM_TAG_TYPE_UINT8) :
      // force value into integer-type (instead of char value)
      out << "i:" << static_cast<uint16_t>(static_cast<uint8_t>(tag_data[index]));
      ++index;
      break;      
    case (Constants::BAM_TAG_TYPE_INT16) :
      out << "i:" << BamTools::UnpackSignedShort(&tag_data[index]);
      index += sizeof(int16_t);
      break;      
    case (Constants::BAM_TAG_TYPE_UINT16) :
      out << "i:" << BamTools::UnpackUnsignedShort(&tag_data[index]);
      index += sizeof(uint16_t);
      break;      
    case (Constants::BAM_TAG_TYPE_INT32) :
      out << "i:" << BamTools::UnpackSignedInt(&tag_data[index]);
      index += sizeof(int32_t);
      break;      
    case (Constants::BAM_TAG_TYPE_UINT32) :
      out << "i:" << BamTools::UnpackUnsignedInt(&tag_data[index]);
      index += sizeof(uint32_t);
      break;      
    case (Constants::BAM_TAG_TYPE_FLOAT) :
      out << "f:" << BamTools::UnpackFloat(&tag_data[index]);
      index += sizeof(float);
      break;      
    case (Constants::BAM_TAG_TYPE_HEX)    : // Fall Through
      ;
    case (Constants::BAM_TAG_TYPE_STRING) :
      out << type << ":";
      while (tag_data[index]) {
        out << tag_data[index];
        ++index;
      }
      ++index;
      break;
    case (Constants::BAM_TAG_TYPE_ARRAY)  :
      char sub_type = tag_data[index];
      out << "B:" << sub_type;
      ++index;
      int array_len = BamTools::UnpackSignedInt(&tag_data[index]);
      index += sizeof(int32_t);
      //      std::cerr << tag_name << ":B:" << sub_type << ":" << array_len << std::endl;
      for(int i=0; i<array_len; i++) {
        switch( sub_type ) {
        case (Constants::BAM_TAG_TYPE_INT8) :
          // force value into integer-type (instead of char value)
          out << "," << static_cast<int16_t>(static_cast<int8_t>(tag_data[index]));
          ++index;
          break;      
        case (Constants::BAM_TAG_TYPE_UINT8) :
          // force value into integer-type (instead of char value)
          out << "," << static_cast<uint16_t>(static_cast<uint8_t>(tag_data[index]));
          ++index;
          break;      
        case (Constants::BAM_TAG_TYPE_INT16) :
          out << "," << BamTools::UnpackSignedShort(&tag_data[index]);
          index += sizeof(int16_t);
          break;      
        case (Constants::BAM_TAG_TYPE_UINT16) :
          out << "," << BamTools::UnpackUnsignedShort(&tag_data[index]);
          index += sizeof(uint16_t);
          break;      
        case (Constants::BAM_TAG_TYPE_INT32) :
          out << "," << BamTools::UnpackSignedInt(&tag_data[index]);
          index += sizeof(int32_t);
          break;      
        case (Constants::BAM_TAG_TYPE_UINT32) :
          out << "," << BamTools::UnpackUnsignedInt(&tag_data[index]);
          index += sizeof(uint32_t);
          break;      
        case (Constants::BAM_TAG_TYPE_FLOAT) :
          out << "," << BamTools::UnpackFloat(&tag_data[index]);
          index += sizeof(float);
          break;
        }              
      }
      break;
    }  
    if ( tag_data[index] == '\0' )
      break;
  }

  out << std::endl;
  return out.str();
}

}
