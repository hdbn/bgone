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
#include "sa2phi.hpp"
#include <iostream>
#include <fstream>

using namespace LZShuffle;

void show2(int *A, int n){
  for(int i = 0; i < n; i++){
    std::cout << std::max(-1, A[i]) << ", ";
  }
  std::cout << std::endl;
}

int main(int argc, char * argv[]){
  std::string s;
  // parse options and read/construct string & suffix array
  unsigned int * sa = (unsigned int*)Init(argc, argv, s);
  int n = s.size(), i;
  unsigned char * us = new unsigned char[n];
  for (i = 0; i < n; i++)us[i] = (unsigned char) s[i];
  s.clear();
  // std::vector<int> phi(n), psv(n);
  // std::cout << "sa[0]=" << sa[0] << std::endl;
  // std::cout << "n=" << n << std::endl;
  unsigned int * phi = new unsigned int[n];

  int first = sa[0];
  phi[sa[0]] = -1;
  for(i = 1; i < n; i++) phi[sa[i]] = sa[i-1];
  // for(i = 0; i < n; i++) if (phi[i] >= 0) std::cout << phi[i] << "["<< s[phi[i]] <<"]"<< std::endl;
  sa2phi(us, sa, n);
  // std::cout << "------- true ---------" << std::endl;
  // show2(phi, n);
  // std::cout << "------- true ---------" << std::endl;
  // show2(sa, n);
  // sa[n/2] = 3;
  for (i = 0; i < n; i++){
    if (i == first) continue;
    if (sa[i] != phi[i]){
      std::cout << "i=" << i << std::endl;
      break;
    }
  }
  delete [] phi;
  delete [] sa;
  delete [] us;
  // std::cout << "Time for phi: " << gettime() - t1 << std::endl;
  // t2 = gettime();
  // for(i = 0; i < n; i++) psv[i] = nsv[i] = -1;
  // for(i = 0; i < n; i++) peakElim(phi[i], i, psv, nsv, -1);
  // std::cout << "Time for pnsv: " << gettime() - t2 << std::endl;

  // ////////////////////////////////////////////////////////////
  // // calculate LZ factorization from text order PSV, NSV
  // ////////////////////////////////////////////////////////////
  // std::vector<std::pair<int,int> > lz;
  // t2 = gettime();
  // lzFromTOPNSV(s, psv, nsv, lz);
  // std::cout << "Time for lz: " << gettime() - t2 << std::endl;
  // std::cout << "# of lz factors: " << lz.size() << std::endl;
  // std::cout << "Total: " << gettime() - t1 << std::endl;
  // if(checkResult){
  //   std::string t = lz2str(lz);
  //   if(s != t) std::cerr << "CHECK: ERROR: mismatch" << std::endl;
  //   else std::cerr << "CHECK: OK" << std::endl;
  // }
  // return 0;
}
