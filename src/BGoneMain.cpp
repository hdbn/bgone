////////////////////////////////////////////////////////////////////////////////
// BGoneMain.cpp
//   lz factorization like one by Karkkainen et.al in CPM 2013.
//   uses only 5n bytes and constant space.
/////////////////////////////////////////////////////////////////////

#include "bgCommon.hpp"
#include "sa2phi.hpp"
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <cassert>
#include <getopt.h>

using namespace LZShuffle;


void usage(int argc, char * argv []){
    std::cout << "Usage  : " << argv[0] << " [options]" << std::endl
	      << "Options: " << std::endl
	      << "  -f iFile : file to process" << std::endl
	      << "  -c       : counts factor only" << std::endl
              << "  -a 1: compute text -> nsv" << std::endl
              << "  -a 2: compute sa   -> nsv" << std::endl
              << "  -a 3: compute text -> nsv (use constant stack)" << std::endl
              << "  -a 4: compute sa   -> nsv (use constant stack)" << std::endl;

}

int main(int argc, char * argv[]){  int ch;
  std::string inFile;
  int algoType = 0;
  bool count_only = false;
  while ((ch = getopt(argc, argv, "a:f:ghc")) != -1) {
    switch (ch) {
    case 'f':
      inFile = optarg;
      break;
    case 'a':
      algoType = atoi(optarg);
      break;
    case 'c':
      count_only = true;
      break;
    case 'g':
      checkResult = true;
      break;

    default:
      usage(argc, argv);
      exit(0);
    }
  }
  if(inFile.empty() || algoType < 1 || algoType > 4){ usage(argc, argv); exit(0); }


  long double t2 = 0;
  
  unsigned int * sa = NULL;
  unsigned int n;
  unsigned char * s;
  long double start_fromS = -1;
  long double start_fromSA = -1;
  std::string InitS;
  std::vector<std::pair<int,int> > * lz;
  std::cerr << "count only: ";
  if (count_only){
    std::cerr << "yes" << std::endl;
    lz = NULL;
  }else{
    std::cerr << "no" << std::endl;
    lz = new std::vector<std::pair<int,int> >();
  }
  if (algoType == 1 || algoType == 3){
    read_text(inFile.c_str(), s, n); // Note size of s is n+1 not n.
    assert(s[n]==0);

    start_fromS = gettime();
    start_fromSA = -1;
    t2 = gettime();
    sa = new unsigned int[n+1];

    ////////////////////////////////////////////////////////////
    // calculate text order PSV, NSV using phi
    ////////////////////////////////////////////////////////////

    if (algoType == 3) text2nsv(s, sa, n+1, CHAR_SIZE, true);
    else text2nsv(s, sa, n+1, CHAR_SIZE, false);
    sa++;
  }else if (algoType == 2 || algoType == 4){
    stringFromFile(inFile, InitS);
    n = InitS.size();

    start_fromS = gettime();
    t2 = gettime();
    sa = (unsigned int *) suffixArray(InitS, (int *)sa, 0);
    std::cerr << "Time for sa: " << gettime() - t2 << std::endl;

    start_fromSA = gettime();
    t2 = gettime();
    s = (unsigned char *) InitS.c_str();
    if (algoType == 4) sa2nsv(s, sa, n, true);
    else sa2nsv(s, sa, n, false);
  }
  std::cerr << "Time for pnsv: " << gettime() - t2 << std::endl;

  ////////////////////////////////////////////////////////////
  // calculate LZ factorization from text order PSV, NSV
  ////////////////////////////////////////////////////////////
  t2 = gettime();
  unsigned int num_factor = lzFromTONSV2(s, sa, n, lz);
  std::cerr << "Time for lz: " << gettime() - t2 << std::endl;
  std::cerr << "Time for total_from_text: " << gettime() - start_fromS << std::endl;
  if (start_fromSA == -1){
  }else{
    std::cerr << "Time for total_from_SA: " << gettime() - start_fromSA << std::endl;
  }
  std::cerr << "# of lz factors: " << num_factor << std::endl;
  if(checkResult){
    std::string t = lz2str(*lz);
    bool ok = true;
    if (n != t.size()) ok = false;
    for (unsigned int i = 0; i < n; i++){
      if(s[i] != (unsigned char)t[i]){
        ok = false;
        break;
      }
    }
    if(!ok) std::cerr << "CHECK: ERROR: mismatch" << std::endl;
    else std::cerr << "CHECK: OK" << std::endl;
  }
  if (lz != NULL) delete lz;
  return 0;
}
