////////////////////////////////////////////////////////////////////////////////
// bgCommon.hpp
////////////////////////////////////////////////////////////////////////////////

#ifndef __LZBG_COMMON_HPP__
#define __LZBG_COMMON_HPP__

#include <string>
#include <vector>
#include <climits>

namespace LZShuffle {
  const unsigned int EMPTY = UINT_MAX; /**< same definition of EMPTY in sa2phi.hpp */

  extern bool checkResult;
  extern bool count_num_factors_only;
  
  enum FLAGS {
    DOUBLE_SA = 1,              /**< allocate double required memory for suffix array */
    SENTINEL = 2                /**< put dummy value to sa[0] and sa[n+1], sa[1..n] is real suffix array of s */
  };

  /** 
   * print usage information
   * 
   * @param argc 
   * @param argv 
   */
  void print_usage(int argc, char * argv []);

  /** 
   * parse options and read/construct string & suffix array
   * allocates memory if sa == 0. return memory for sa.
   * 
   * @param argc 
   * @param argv 
   * @param s input text
   * @param f option for suffix array
   * @param sa pointer of suffix array
   * 
   * @return sa // does it need?
   */
  int * Init(int argc, char * argv[], std::string & s, 
	     unsigned int f = 0, int *sa = 0);

  /** 
   * read file fileName into string s.
   * 
   * @param fileName file name
   * @param s string that input text will be stored
   */
  void stringFromFile(const std::string & fileName, std::string & s);

  /** 
   * calculate suffix array sa from s. allocates memory if sa == 0
   * 
   * @param s input text
   * @param sa suffix array
   * @param f option for suffix array
   * 
   * @return sa
   */
  int * suffixArray(const std::string & s, int * sa, unsigned int f);

  /** 
   * read suffix array of string s from fname + '.sa' into sa.
   * if safname does not exist or seems invalid, create suffix array and save it
   * to fname + '.sa'
   * allocates memory if *sa == 0
   * 
   * @param s input text
   * @param fname file name
   * @param sa suffix array
   * @param f option for suffix array
   * 
   * @return 
   */
  int * saFromFile(const std::string & s, const std::string & fname, 
		   int * sa, unsigned int f);


  ////////////////////////////////////////////////////////////////////////////////
  // computing lz from PSV, NSV in text order
  ////////////////////////////////////////////////////////////////////////////////

  /** 
   * compute LZ factorization by using PSV and NSV
   * 
   * @param s input text
   * @param psv psv array
   * @param nsv nsv array
   * @param lz lz factors
   * 
   * @return number of LZ factors
   */
  unsigned int lzFromTOPNSV(const std::string & s, 
		    const int * psv,
		    const int * nsv,
		    std::vector<std::pair<int,int> > * lz);


  /** 
   * compute LZ factorization by using NSV only
   * the algorithm is based on the paper by Karkkainen et.al in CPM 2013.
   * 
   * @param s input text
   * @param nsv nsv array
   * @param n the length of s
   * @param lz lz factors
   * 
   * @return number of LZ factors
   */
  unsigned int lzFromTONSV2(const unsigned char * s, 
                            unsigned int * nsv,
                            unsigned int n,
                            std::vector<std::pair<int,int> > * lz);


  /** 
   * recover string from lz factorization
   * 
   * @param lz lz factors
   * 
   * @return string represented by lz factors.
   */
  std::string lz2str(const std::vector<std::pair<int,int> > & lz);

  /** 
   * for experiment
   * 
   * @return time
   */
  double gettime();

  // read text by Simon Puglisi
  // Copyright (c) 2013 Juha Karkkainen, Dominik Kempa and Simon J. Puglisi
  void read_text(const char *filename, unsigned char* &text, unsigned int &length);
};
#endif//__LZBG_COMMON_HPP__
