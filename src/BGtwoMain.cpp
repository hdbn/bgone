////////////////////////////////////////////////////////////////////////////////
// bgtMain.cpp
//   lz factorization via PSV_text and NSV_text using peak elimination
//   uses 13N bytes space
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

#include "bgCommon.hpp"
#include <iostream>
#include <fstream>

using namespace LZShuffle;

int main(int argc, char * argv[]){

  std::string s;

  // parse options and read/construct string & suffix array
  int * sa = Init(argc, argv, s);
  std::vector<std::pair<int,int> > * lz;
  std::cerr << "count only: ";
  if (count_num_factors_only){
    std::cerr << "yes" << std::endl;
    lz = NULL;
  }else{
    std::cerr << "no" << std::endl;
    lz = new std::vector<std::pair<int, int> >();
  }
  long double time_sa_start = gettime();
  // std::cout << "Time for sa: " << gettime() - t1 << std::endl;
  
  ////////////////////////////////////////////////////////////
  // calculate text order PSV, NSV using phi
  ////////////////////////////////////////////////////////////
  long double time_pnsv_start = gettime();
   int n = s.size(), i;
  int sa_last = sa[n-1];
  // std::vector<int> phi(n), psv(n);
  int * phi = new int[n];
  phi[sa[0]] = -1;

  for(i = 1; i < n; i++) {phi[sa[i]] = sa[i-1];}
  int prev = -1;
  int cur = sa_last;
  while (true){
    const int next = phi[cur];
    while (cur < prev){
      sa[prev] = cur; // sa[prev] == psv[prev]
      prev = phi[prev]; // nsv[prev]
    }
    if (cur < 0) break;
    phi[cur] = prev; // phi[cur] == nsv[cur]
    prev = cur;
    cur = next;
  }
  // Now, phi is nsv, sa is psv array

  // for(i = 0; i < n; i++) psv[i] = nsv[i] = -1;
  std::cerr << "Time for pnsv: " << gettime() - time_pnsv_start << std::endl;

  ////////////////////////////////////////////////////////////
  // calculate LZ factorization from text order PSV, NSV
  ////////////////////////////////////////////////////////////
  long double time_lz_start = gettime();
  unsigned int nfactor = lzFromTOPNSV(s, sa, phi, lz);
  std::cerr << "Time for lz: " << gettime() - time_lz_start << std::endl;
  // std::cout << "# of lz factors: " << lz.size() << std::endl;
  std::cerr << "Time for total_from_SA: " << gettime() - time_sa_start << std::endl;
  std::cerr << "# of lz factors: " << nfactor << std::endl;
  // std::cerr << "Total: " << gettime() - t1 << std::endl;
  if(checkResult){
    std::string t = lz2str(*lz);
    if(s != t) std::cerr << "CHECK: ERROR: mismatch" << std::endl;
    else std::cerr << "CHECK: OK" << std::endl;
  }
  return 0;
}
