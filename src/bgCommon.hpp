////////////////////////////////////////////////////////////////////////////////
// bgCommon.hpp
//   common routines for lz factorization BGS, iBGS, BGL, iBGL, BGT, iBGT
////////////////////////////////////////////////////////////////////////////////
// Copyright 2012 Hideo Bannai & Keisuke Goto
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
////////////////////////////////////////////////////////////////////////////////

#ifndef __LZBG_COMMON_HPP__
#define __LZBG_COMMON_HPP__

#include <string>
#include <vector>
#include <climits>

namespace LZShuffle {
  const unsigned int EMPTY = UINT_MAX; // same definition of EMPTY in sa2phi.hpp

  extern bool checkResult;
  extern bool count_num_factors_only;
  
  enum FLAGS {
    DOUBLE_SA = 1, // allocate double required memory for suffix array
    SENTINEL = 2 // put dummy value to sa[0] and sa[n+1], sa[1..n] is real suffix array of s
  };
  
  // print usage information
  void print_usage(int argc, char * argv []);

  // parse options and read/construct string & suffix array
  // allocates memory if sa == 0. return memory for sa.
  int * Init(int argc, char * argv[], std::string & s, 
	     unsigned int f = 0, int *sa = 0);

  // read file fileName into string s.
  void stringFromFile(const std::string & fileName, std::string & s);

  // calculate suffix array sa from s. allocates memory if sa == 0
  // return *sa
  int * suffixArray(const std::string & s, int * sa, unsigned int f);

  // read suffix array of string s from fname + '.sa' into sa.
  // if safname does not exist or seems invalid, create suffix array and save it
  // to fname + '.sa'
  // allocates memory if *sa == 0
  // return *sa
  int * saFromFile(const std::string & s, const std::string & fname, 
		   int * sa, unsigned int f);


  ////////////////////////////////////////////////////////////////////////////////
  // computing lz from PSV, NSV in text order
  ////////////////////////////////////////////////////////////////////////////////

  // for bgt9
  unsigned int lzFromTOPNSV(const std::string & s, 
		    const int * psv,
		    const int * nsv,
		    std::vector<std::pair<int,int> > * lz);

  unsigned int lzFromTONSV2(const unsigned char * s, 
                            unsigned int * nsv,
                            unsigned int n,
                            std::vector<std::pair<int,int> > * lz);


  // recover string from lz factorization
  std::string lz2str(const std::vector<std::pair<int,int> > & lz);

  double gettime();

  // read text by Simon Puglisi
  void read_text(const char *filename, unsigned char* &text, unsigned int &length);
};
#endif//__LZBG_COMMON_HPP__
