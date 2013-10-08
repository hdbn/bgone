////////////////////////////////////////////////////////////////////////////////
// sa2phi_test.cpp
//   test file for sa2phi.cpp
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
}
