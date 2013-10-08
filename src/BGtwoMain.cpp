////////////////////////////////////////////////////////////////////////////////
// BGtwoMain.cpp
//   lz factorization via PSV_text and NSV_text using peak elimination
//   uses 9N bytes space
///////////////////////////////////////////////////////////////////////

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
  
  ////////////////////////////////////////////////////////////
  // calculate text order PSV, NSV using phi
  ////////////////////////////////////////////////////////////
  long double time_pnsv_start = gettime();
   int n = s.size(), i;
  int sa_last = sa[n-1];
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

  std::cerr << "Time for pnsv: " << gettime() - time_pnsv_start << std::endl;

  ////////////////////////////////////////////////////////////
  // calculate LZ factorization from text order PSV, NSV
  ////////////////////////////////////////////////////////////
  long double time_lz_start = gettime();
  unsigned int nfactor = lzFromTOPNSV(s, sa, phi, lz);
  std::cerr << "Time for lz: " << gettime() - time_lz_start << std::endl;
  std::cerr << "Time for total_from_SA: " << gettime() - time_sa_start << std::endl;
  std::cerr << "# of lz factors: " << nfactor << std::endl;
  if(checkResult){
    std::string t = lz2str(*lz);
    if(s != t) std::cerr << "CHECK: ERROR: mismatch" << std::endl;
    else std::cerr << "CHECK: OK" << std::endl;
  }
  return 0;
}
