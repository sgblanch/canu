#include <vector>
#include <stdint.h>
#include <cstring>
#include <string>
#include <algorithm>
#include "Alignment.H"
#include "SimpleAligner.H"
#include "assert.h"

#include "ssw_cpp.H"
#include <sstream>

SimpleAligner::SimpleAligner() {
}

void SimpleAligner::align(dagcon::Alignment &aln, double errorRate, bool useND) {
   if (useND) alignND(aln, errorRate);
   else alignSW(aln, errorRate);
}

void SimpleAligner::alignND(dagcon::Alignment &aln, double errorRate) {
  NDalignment::NDalignResult ndaln;
  int tolerance = std::min(200, (int)round(0.5 * errorRate * (aln.qstr.size() + aln.tstr.size())));
//fprintf(stderr, "Aligning at %f strings of %d %d tolerance must be %d\n", errorRate, aln.qstr.size(), aln.tstr.size(), tolerance);
  bool aligned = NDalignment::align(aln.qstr.c_str(), aln.qstr.size(), aln.tstr.c_str(), aln.tstr.size(), tolerance, true, ndaln);

//fprintf(stderr, "Alignment is %d and strings are %d, %d\n",((int)aligned), ndaln._dist, ndaln._size);
  if (((double) ndaln._dist / (double) ndaln._size) > errorRate) {
     aligned = false;
  }

  if (aligned) {
     aln.start += ndaln._tgt_bgn;
     aln.end = aln.start + ndaln._tgt_end;
     aln.start++;
     aln.qstr = std::string(ndaln._qry_aln_str);
     aln.tstr = std::string(ndaln._tgt_aln_str);
  } else {
     aln.start = aln.end = 0;
     aln.qstr = std::string();
     aln.tstr = std::string();
  }

  if (aln.qstr.length() != aln.tstr.length()) {
fprintf(stderr, "Found a bad alignment with inconsistent lengths %d %d\n", aln.qstr.length(), aln.tstr.length());
     aln.start = aln.end = 0;
     aln.qstr = std::string();
     aln.tstr = std::string();
  }


//fprintf(stderr, "Aligned sequences of length %d to %d with error rate %f and found an alignment of %d errors (%d length). The positions %d - %d and %d - %d\n", aln.tstr.size(), aln.qstr.size(), errorRate, ndaln._dist, ndaln._size, ndaln._tgt_bgn, ndaln._tgt_end, ndaln._qry_bgn, ndaln._qry_end);
  assert(aln.qstr.length() == aln.tstr.length());
}

void SimpleAligner::alignSW(dagcon::Alignment &aln, double errorRate) {
  StripedSmithWaterman::Aligner aligner(1, 2, 2, 1);
  StripedSmithWaterman::Filter filter;
  StripedSmithWaterman::Alignment alignment;
  bool aligned = aligner.Align(aln.qstr.c_str(), aln.tstr.c_str(), aln.tstr.size(), filter, &alignment);

  std::istringstream parser(alignment.cigar_string);
  std::string tgtStr;
  std::string qryStr;

  char c;
  uint32_t len;
  int posT = alignment.ref_begin;
  int posQ = alignment.query_begin;
  int total = 0;
  int matches = 0;
//fprintf(stderr, "Aligned sequences with positions %d - %d and %d - %d cigar %s\n", alignment.ref_begin, alignment.ref_end, alignment.query_begin, alignment.query_end, alignment.cigar_string.c_str());
  while (parser >> len >> c) {
     total += len;
     switch ( c) {
        case 'M':
        case '=':
	   matches++;
        case 'X':
           tgtStr.append(aln.tstr.substr(posT, len));
           qryStr.append(aln.qstr.substr(posQ, len));
           posT+=len;
           posQ+=len;
           break;
        case 'I':
           tgtStr.append(std::string(len, '-'));
           qryStr.append(aln.qstr.substr(posQ, len));
           posQ+=len;
           break;
        case 'D':
           tgtStr.append(aln.tstr.substr(posT, len));
           qryStr.append(std::string(len, '-'));
           posT+=len;
           break;
        case 'S':
           break;
        default:
           fprintf(stderr, "Unknown character %c\n", c);
     }
  }
  if ((double)alignment.mismatches / (double) total > (errorRate)) {
     aligned = false;
  }

//fprintf(stderr, "Aligned sequences of length %d to %d with error rate %f and found an alignment of %d errors (%d length). The positions %d - %d and %d - %d\n", aln.tstr.size(), aln.qstr.size(), errorRate, alignment.mismatches, total, alignment.ref_begin, alignment.ref_end, alignment.query_begin, alignment.query_end);

  if (aligned) {
     aln.start += alignment.ref_begin;
     aln.end = aln.start + alignment.ref_end;
     aln.start++;

//     fprintf(stderr, "I have all done with length %d - %d and originally they were %d - %d total len %d\n", tgtStr.length(), qryStr.length(), aln.tstr.length(), aln.qstr.length(), total);
     aln.qstr = qryStr;
     aln.tstr = tgtStr;
  } else {
     aln.start = aln.end = 0;
     aln.qstr = std::string();
     aln.tstr = std::string();
  }

  assert(aln.qstr.length() == aln.tstr.length());
}
