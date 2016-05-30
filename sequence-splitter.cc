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
#include <api/BamConstants.h>
#include <api/BamReader.h>
#include "sequence-splitter.h"

namespace Constants = ::BamTools::Constants;
using BamTools::BamAlignment;
using BamTools::CigarOp;
using std::string;
using std::vector;


namespace fmu_tools {

// print BamAlignment in SAM format

string SequenceSplitter::ConvertToSam(const BamAlignment& d) {
  // tab-delimited
  // <QNAME> <FLAG> <RNAME> <POS> <MAPQ> <CIGAR> <MRNM> <MPOS> <ISIZE> <SEQ> <QUAL> [ <TAG>:<VTYPE>:<VALUE> [...] ]
  // write name & alignment flag
  std::ostringstream m_out;
  m_out << d.Name << "\t" << d.AlignmentFlag << "\t";

  /*
  // write reference name
  if ( (a.RefID >= 0) && (a.RefID < (int)m_references.size()) ) 
    m_out << m_references[a.RefID].RefName << "\t";
  else 
  */
  m_out << "*\t";
  
  // write position & map quality
  m_out << d.Position+1 << "\t" << d.MapQuality << "\t";
    
  // write CIGAR
  const vector<CigarOp>& cigarData = d.CigarData;
  if ( cigarData.empty() ) 
    m_out << "*\t";
  else {
    vector<CigarOp>::const_iterator cigarIter = cigarData.begin();
    vector<CigarOp>::const_iterator cigarEnd  = cigarData.end();
    for ( ; cigarIter != cigarEnd; ++cigarIter ) {
      const CigarOp& op = (*cigarIter);
      m_out << op.Length << op.Type;
    }
    m_out << "\t";
  }
  /*
  // write mate reference name, mate position, & insert size
  if ( d.IsPaired() && (d.MateRefID >= 0) && (d.MateRefID < (int)m_references.size()) ) {
        if ( a.MateRefID == a.RefID )
            m_out << "=\t";
        else
           m_out << m_references[a.MateRefID].RefName << "\t";
        m_out << a.MatePosition+1 << "\t" << a.InsertSize << "\t";
    } 
    else
  */
  m_out << "*\t0\t0\t";
    
  // write sequence
  if ( d.QueryBases.empty() )
    m_out << "*\t";
  else
    m_out << d.QueryBases << "\t";
    
  // write qualities
  if ( d.Qualities.empty() || (d.Qualities.at(0) == (char)0xFF) )
    m_out << "*";
  else
    m_out << d.Qualities;
    
  // write tag data
  const char* tagData = d.TagData.c_str();
  const size_t tagDataLength = d.TagData.length();
    
  size_t index = 0;
  while ( index < tagDataLength ) {
    // write tag name   
    string tagName = d.TagData.substr(index, 2);
    m_out << "\t" << tagName << ":";
    index += 2;        
    // get data type
    char type = d.TagData.at(index);
    ++index;
    switch ( type ) {
    case (Constants::BAM_TAG_TYPE_ASCII) :
      m_out << "A:" << tagData[index];
      ++index;
      break;
      
    case (Constants::BAM_TAG_TYPE_INT8) :
      // force value into integer-type (instead of char value)
      m_out << "i:" << static_cast<int16_t>(tagData[index]);
      ++index;
      break;
      
    case (Constants::BAM_TAG_TYPE_UINT8) :
      // force value into integer-type (instead of char value)
      m_out << "i:" << static_cast<uint16_t>(tagData[index]);
      ++index;
      break;
      
    case (Constants::BAM_TAG_TYPE_INT16) :
      m_out << "i:" << BamTools::UnpackSignedShort(&tagData[index]);
      index += sizeof(int16_t);
      break;
      
    case (Constants::BAM_TAG_TYPE_UINT16) :
      m_out << "i:" << BamTools::UnpackUnsignedShort(&tagData[index]);
      index += sizeof(uint16_t);
      break;
      
    case (Constants::BAM_TAG_TYPE_INT32) :
      m_out << "i:" << BamTools::UnpackSignedInt(&tagData[index]);
      index += sizeof(int32_t);
      break;
      
    case (Constants::BAM_TAG_TYPE_UINT32) :
      m_out << "i:" << BamTools::UnpackUnsignedInt(&tagData[index]);
      index += sizeof(uint32_t);
      break;
      
    case (Constants::BAM_TAG_TYPE_FLOAT) :
      m_out << "f:" << BamTools::UnpackFloat(&tagData[index]);
      index += sizeof(float);
      break;
      
    case (Constants::BAM_TAG_TYPE_HEX)    : // fall-through
    case (Constants::BAM_TAG_TYPE_STRING) :
      m_out << type << ":";
    while (tagData[index]) {
      m_out << tagData[index];
      ++index;
    }
    ++index;
    break;
    }
    
    if ( tagData[index] == '\0' )
      break;
  }
  m_out << std::endl;
  return m_out.str();
}


} // fmu_tools
