////////////////////////////////////////////////////////////////////////////////
// bgCommon.cpp
////////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#include <sys/time.h>
#include <cassert>
#include <cstdlib>
#include "bgCommon.hpp"
#include "divsufsort.h"

#define PSV(i) pnsv[(i << 1)]
#define NSV(i) pnsv[(i << 1)+1]

#define LPS(i) lpspo[(i << 1)]
#define PREVOCC(i) lpspo[(i << 1)+1]

namespace LZShuffle {

  bool checkResult = false;
  bool count_num_factors_only = false;

  ////////////////////////////////////////////////////////////
  // parse options and read/construct string & suffix array
  ////////////////////////////////////////////////////////////
  
  int * Init(int argc, char * argv[], std::string & s, unsigned int f, int * sa){
    int ch;
    std::string inFile, saFile;
    bool useSAcache = false;
    optind = 1;
    optreset = 1; // for using getopt at multiple times.
    while ((ch = getopt(argc, argv, "f:xghc")) != -1) {
      switch (ch) {
      case 'f':
	inFile = optarg;
	break;
      case 'x':
	useSAcache = true;
	break;
      case 'g':
	checkResult = true;
	break;
      case 'c':
	count_num_factors_only = true;
	break;
      default:
	print_usage(argc, argv);
	exit(0);
      }
    }
    if(inFile.empty()){ print_usage(argc, argv); exit(0); }
    
    argc -= optind;
    argv += optind;
    ////////////////////////////////////////////////////////////
    
    stringFromFile(inFile, s);
    if(useSAcache){
      sa = saFromFile(s, inFile, sa, f);
    } else {
      sa = suffixArray(s, sa, f);
    }
    return(sa);
  }

  void print_usage(int argc, char * argv []){
    std::cout << "Usage  : " << argv[0] << " [options]" << std::endl
	      << "Options: " << std::endl
	      << "  -f iFile : file to process" << std::endl
	      << "  -x       : use iFile + '.sa' for suffix array cache" << std::endl;
    return;
  }

  ////////////////////////////////////////////////////////////
  //read input string from file
  ////////////////////////////////////////////////////////////
  void stringFromFile(const std::string & fileName, std::string & s){
    struct stat st;
    size_t fileSize;
    if(stat(fileName.c_str(), &st)){
      std::cerr << "failed to stat file: " << fileName << std::endl;
      return;
    }
    std::ifstream ifs(fileName.c_str(), std::ios::in | std::ios::binary);
    if(!ifs){
      std::cerr << "failed to read file: " << fileName << std::endl;
      return;
    }
    fileSize = st.st_size;
    if(fileSize != static_cast<size_t>(static_cast<int>(fileSize))){
      std::cerr << "ERROR: The file size is too big to fit in int. Cannot process." << std::endl;
      return;
    }
    s.resize(fileSize);
    ifs.read(reinterpret_cast<char*>(&s[0]), fileSize);
  }


  int * suffixArray(const std::string & s, int * sa, unsigned int f){
    int plus_size = 0;
    if (f & SENTINEL) plus_size = 2;
    if(sa == 0){
      size_t sasize = (f & DOUBLE_SA) ? s.size() * 2 : s.size();
      sa = new int[sasize+plus_size];
    }
    std::cerr << "Building suffix array..." << std::flush;
    // std::cout << (f == SENTINEL) << std::endl;
    divsufsort(reinterpret_cast<const unsigned char *>(s.c_str()), sa+(f == SENTINEL), s.size());
    std::cerr << "done" << std::endl;
    return (sa);
  }

  int * saFromFile(const std::string & s, const std::string & fname, int * sa, unsigned int f){
    ////////////////////////////////////////////////////////////
    // check if suffix array file if it exists
    ////////////////////////////////////////////////////////////
    bool remake = false;
    struct stat st1, st2;
    std::string safname = fname + ".sa";
    if(stat(fname.c_str(), &st1)) remake = true;
    if(stat(safname.c_str(), &st2)) remake = true;
    if(st1.st_mtime >= st2.st_mtime){
      remake = true;
    }
    std::ifstream sfs(safname.c_str(), std::ios::in | std::ios::binary);
    int plus_size = 0;
    if (f & SENTINEL) plus_size = 2;
    if(!sfs){
      remake = true;
    } else if(sfs){
      const size_t fileSize = st2.st_size;
      if (fileSize != sizeof(int) * s.size()){
	remake = true;
      } else {
	std::cerr << "reading suffix array from: " << safname << std::flush;
	if(sa == 0){
	  size_t sasize = (f & DOUBLE_SA) ? (fileSize / sizeof(int)) * 2 : fileSize / sizeof(int);
	  sa = new int[sasize+plus_size];
	}
	sfs.read(reinterpret_cast<char*>(sa+((f&SENTINEL)>0)), fileSize);
	std::cerr << " ...done" << std::endl;
      }
    }
    ////////////////////////////////////////////////////////////
    // construct suffix array and cache it to file
    ////////////////////////////////////////////////////////////
    if(remake){
      std::cerr << "suffix array file: " << safname << " not found, invalid or out of date." << std::endl;
      sa = suffixArray(s, sa, f);
      std::cerr << "Saving suffix array to file..." << std::flush;
      std::ofstream ofs(safname.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
      ofs.write(reinterpret_cast<const char*>(sa+((f&SENTINEL)>0)), sizeof(int) * s.size());
      std::cerr << "done" << std::endl;
    }
    return(sa);
  }

  inline int naiveLCP(const unsigned char * const x, int i, int j, int n){
    int l = 0;
    // std::cout << "(i, j, n)=(" << i << ", " << j << ", " << n << ")" << std::endl;
    assert(0 <= i && i < j);
    // if (j < i) std::swap(i, j);
    while(j < n && x[i++] == x[j++]){ l++; }
    return l;
  }

  std::pair<int, int> factorTOPNSV2(const unsigned char * s,
                                   unsigned int n,
                                   unsigned int p,
                                   unsigned int psvi,
                                   unsigned int nsvi){
    int plen = (psvi > p) ? 0 : naiveLCP(s, psvi, p, n);
    int nlen = (nsvi > p) ? 0 : naiveLCP(s, nsvi, p, n);
    // int plen = (psvi >= n) ? 0 : naiveLCP(s, psvi, p, n);
    // int nlen = (nsvi >= n) ? 0 : naiveLCP(s, nsvi, p, n);
    // std::cout << std::endl << s << std::endl;
    // std::cout << "(plen, nlen)=(" << plen << "," << nlen << ") s[" << p << "]=" << s[p] << std::endl;
    if (plen > nlen) return std::make_pair(plen, psvi);
    if (nlen > 0) return std::make_pair(nlen ,nsvi);
    else return std::make_pair(0, s[p]);
  }
  std::pair<int, int> factorTOPNSVOpt2(const unsigned char * s,
                                       unsigned int n,
                                       int p,
                                       unsigned int psvi,
                                       unsigned int nsvi){
    int prevPos;
    int common_len = 0;
    if (psvi >= n && nsvi >= n){
      common_len = naiveLCP(s, psvi, p, n);
      int plen = naiveLCP(s, psvi + common_len, p + common_len, n);
      if (plen > 0){
        common_len += plen;
        prevPos = psvi;
      }else{
        common_len += naiveLCP(s, nsvi + common_len, p + common_len, n);
        prevPos = nsvi;
      }
    }else if (psvi < n || nsvi < n){
      prevPos = psvi < n ? psvi : nsvi;
      common_len = naiveLCP(s, prevPos, p, n);
    }
    if(common_len > 0){
      return std::make_pair(common_len, prevPos);
    } else {
      return std::make_pair(0, s[p]);
    }
  }

  unsigned int lzFromTOPNSV(const std::string & s, 
		    const int * psv,
		    const int * nsv,
		    std::vector<std::pair<int,int> > * lz){
    
    ////////////////////////////////////////////////////////////
    // main LZ factorization
    ////////////////////////////////////////////////////////////
    size_t p = 1;
    if (lz) lz->push_back(std::make_pair(0, s[0]));
    unsigned int num_factor = 1;
    while(p < s.size()){
      int prevPos = psv[p];
      int lpf = (psv[p] < 0) ? 0 : naiveLCP((unsigned char *)s.c_str(), psv[p], p, s.size());
      int nlen = (nsv[p] < 0) ? 0 : naiveLCP((unsigned char *)s.c_str(), nsv[p], p, s.size());
      if(nlen > lpf){ lpf = nlen; prevPos = nsv[p]; }
      if(lpf > 0){
	if (lz) lz->push_back(std::make_pair(lpf,prevPos));
	p += lpf;
      } else {
	if (lz) lz->push_back(std::make_pair(0, s[p]));
	p++;
      }
      num_factor++;
    }
    return num_factor;
    ////////////////////////////////////////////////////////////
  }

  void lzFromTOPNSV(const std::string & s, 
		    const int * pnsv,
		    std::vector<std::pair<int,int> > & lz){
    
    ////////////////////////////////////////////////////////////
    // main LZ factorization
    ////////////////////////////////////////////////////////////
    size_t p = 1;
    lz.clear();
    lz.push_back(std::make_pair(0, s[0]));
    while(p < s.size()){
      int prevPos = PSV(p);
      int lpf = (PSV(p) < 0) ? 0 : naiveLCP((unsigned char *)s.c_str(), PSV(p), p, s.size());
      int nlen = (NSV(p) < 0) ? 0 : naiveLCP((unsigned char *)s.c_str(), NSV(p), p, s.size());
      if(nlen > lpf){ lpf = nlen; prevPos = NSV(p); }
      if(lpf > 0){
	lz.push_back(std::make_pair(lpf,prevPos));
	p += lpf;
      } else {
	lz.push_back(std::make_pair(0, s[p]));
	p++;
      }
    }
    return;
    ////////////////////////////////////////////////////////////
  }

  // void lzFromTONSV2(const unsigned char * s, 
  //                   unsigned int * nsv,
  //                   unsigned int n,
  //                   std::vector<std::pair<int,int> > & lz){
  unsigned int lzFromTONSV2(const unsigned char * s, 
                    unsigned int * nsv,
                    unsigned int n,
                    std::vector<std::pair<int,int> > * lz){
    
    ////////////////////////////////////////////////////////////
    // main LZ factorization
    ////////////////////////////////////////////////////////////
    unsigned int p = 1;
    // lz->clear();
    if (lz != NULL) lz->push_back(std::make_pair(0, s[0]));
    // int l, prevPos;
    unsigned int nsvi, psvi;
    unsigned int next = 1;
    // std::cout << "hoge" << std::endl;
    unsigned int imax = 0;
    unsigned int num_factor = 1;
    // unsigned int EMPTY = UINT_MAX;
    nsv[0] = EMPTY;
    // std::cout << nsv[0] <<", "<< nsv[1] <<", "<< nsv[2] << std::endl;
    while(p < n){
      // std::cout << "p=" << p <<", nsv[p]="<< nsv[p] <<std::endl;
      assert(nsv[p] == EMPTY || nsv[p] < n);
      nsvi = nsv[p];
      if (nsvi != EMPTY){
        psvi = nsv[nsvi];
        nsv[nsvi] = p;
      }else{
        psvi = imax;
        imax = p;
        // nsvi = EMPTY; // because nsvi is EMPTY
      }
      if (p == next){
        // std::cout << "(p=" << p << ", psvi="  << psvi << ", nsvi="<< nsvi << ") ";
        const std::pair<int, int> f = factorTOPNSV2(s, n, p, psvi, nsvi);
        // std::cout << "(" << f.first << "," << f.second << ") ";
        if (lz != NULL) lz->push_back(f);
        next = p+std::max(1, f.first);
        num_factor++;
      }
      nsv[p] = psvi;
      // nsv[nsvi] = p;
      p++;
    }
    return num_factor;
    ////////////////////////////////////////////////////////////
  }
  // void lzFromTOPSV(const std::string & s, 
  //       	    int * psv,
  //       	    std::vector<std::pair<int,int> > & lz){
    
  //   ////////////////////////////////////////////////////////////
  //   // main LZ factorization
  //   ////////////////////////////////////////////////////////////
  //   size_t p = 1;
  //   lz.clear();
  //   lz.push_back(std::make_pair(0, s[0]));
  //   // int l, prevPos;
  //   int nsvi, psvi;
  //   size_t next = 1;
  //   int imax = 0;
  //   while(p < s.size()){
  //     assert(psv[p] < (int)s.size());
  //     psvi = psv[p];
  //     if (psvi > 0){
  //       nsvi = psv[psvi];
  //       psv[psvi] = p;
  //     }else{
  //       // nsvi = psv[s.size()];
  //       nsvi = imax;
  //       // psv[imax] = p;
  //       imax = p;
  //       // psvi = s.size();
  //     }
  //     if (p == next){
  //       const std::pair<int, int> f = factorTOPNSVOpt(s, p, psvi, nsvi);
  //       // std::cout << "(" << f.first << "," << f.second << ") ";
  //       lz.push_back(f);
  //       next = p+std::max(1, f.first);
  //     }
  //     psv[p] = nsvi;
  //     // nsv[nsvi] = p;
  //     p++;
  //   }
  //   return;
  //   ////////////////////////////////////////////////////////////
  // }


  std::string lz2str(const std::vector<std::pair<int,int> > & lz){
    std::string s;
    size_t i;
    int j;
    for(i = 0; i < lz.size(); i++){
      if(lz[i].first == 0){
	s.push_back(static_cast<char>(lz[i].second));
      } else {
	for(j = 0; j < lz[i].first; j++){
	  s.push_back(s[lz[i].second+j]);
	}
      }
    }
    return(s);
  }

  double gettime(){
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return (((double)tv.tv_sec)
	    + ((double)tv.tv_usec)*1e-6);
  }

  void read_text(const char *filename, unsigned char* &text, unsigned int &length) {
    std::fstream f(filename, std::fstream::in);
    if (f.fail()) {
      std::cerr << "\nError: cannot open file " << filename << "\n";
      std::exit(EXIT_FAILURE);
    }

    f.seekg(0, std::ios_base::end);
    length = f.tellg();
    f.seekg (0, std::ios_base::beg);

    text = new unsigned char[length+1];
    if (!text) {
      std::cerr << "\nError: allocation of " << length << " bytes failed\n";
      std::exit(EXIT_FAILURE);
    }

    std::cerr << "Reading the file " << filename << " (" << length << " bytes)... ";
    f.read((char *)text, length);
    if (!f) {
      std::cerr << "\nError: failed to read " << length << " bytes from file "
                << filename << ". Only " << f.gcount() << " could be read\n";
      std::exit(EXIT_FAILURE);
    }
    std::cerr << std::endl;
    f.close();
  }
}

