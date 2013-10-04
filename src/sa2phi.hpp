////////////////////////////////////////////////////////////////////////////////
// sa2phi.hpp
////////////////////////////////////////////////////////////////////////////////
#pragma once
#include <string>
#include <climits>

const unsigned int EMPTY = UINT_MAX;
const int CHAR_SIZE = 256;
const bool LTYPE = false;
const bool STYPE = true;

/** 
 * rewrite from SA to Left Most S-type(LMS) suffix array
 * we jsut do sampling LMS-suffixes from SA
 * 
 * @param s input text
 * @param sa suffix array of s
 * @param n length of s
 * @param n1 number of Left Most S-type(LMS) suffixes
 */
void sa2LMSsa(const unsigned char * s, unsigned int * sa, unsigned int n, unsigned int & n1);

/** 
 * rewrite from SA to successor link of SA in-place
 * 
 * @param s input text
 * @param sa suffix array of s
 * @param n length of s
 * @param sa0 position that lexicographically smallest suffix starts
 */
void sa2phi(const unsigned char * s, unsigned int * sa, unsigned int n);

/** 
 * rewrite from SA to successor link of LMS-SA in-place
 * 
 * @param s input text
 * @param sa suffix array of s
 * @param n length of s
 * @param sa0 position that lexicographically smallest suffix starts
 */
void sa2LMSphi(const unsigned char * s, unsigned int * sa, unsigned int n, int sa0);

/** 
 * induced sort for L-type suffixes on SA from left to right
 * 
 * @param s input text
 * @param phi phi array
 * @param n length of s
 * @param sa0 a position that lexicographically smallest suffix starts
 * @param Lbkts suffix pointers that L-intervals start
 * @param Lbkte suffix pointers that L-intervals end
 */
void inducePhiL(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                std::vector<unsigned int> & Lbkts, std::vector<unsigned int> & Lbkte);

/** 
 * induced sort for S-type suffixes on SA from right to left
 * 
 * @param s input text
 * @param phi phi array
 * @param n length of s
 * @param sa0 position that lexicographically smallest suffix starts
 * @param Lbkts suffix pointers that L-intervals start
 * @param Lbkte suffix pointers that L-intervals end
 */
void inducePhiS(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte);

/** 
 * rewrite from nsv to phi array in-place
 * 
 * @param phi nsv array
 * @param last_idx position that lexicographically largest suffix starts
 */
void nsvTOFromPhi(unsigned int * phi, unsigned int last_idx);


/** 
 * Peak Elimination
 * XXXXXXXXX UNUSED FUNCTION XXXXXXXX
 * 
 * @param j position of suffix
 * @param i position of suffix
 * @param nsv nsv array
 * @param bot 
 */
void peakElim_ow(int j, int i, unsigned int * nsv, int bot);

/** 
 * rewrite from L-type phi array to nsv array by using a stack of constant size
 * 
 * @param s input text
 * @param phi L-type phi array
 * @param n length of s
 * @param sa0 unused variable?
 * @param Lbkts suffix pointers that L-intervals start
 * @param Lbkte suffix pointers that L-intervals end
 */
void phiL2nsv_CS(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                 const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte);


/** 
 * rewrite from L-type phi array to nsv array
 * 
 * @param s input text
 * @param phi L-type phi array
 * @param n length of s
 * @param sa0 unused variable?
 * @param Lbkts suffix pointers that L-intervals start
 * @param Lbkte suffix pointers that L-intervals end
 */
void phiL2nsv(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
              const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte);

/** 
 * Peak Elimination
 * we only compare suffix j and lexicographically larger suffixes than suffix j
 * 
 * @param j position of suffix
 * @param i position of suffix
 * @param nsv nsv[1..j] is phi array, nsv[j+1..n] is nsv array
 * 
 * @return nsv(j) if nsv(j) exists, EMPTY otherwise
 */
int peakElimRight(unsigned int j, unsigned int i, const unsigned int * nsv);

/** 
 * rewrite SA to nsv array
 * 
 * @param s input text
 * @param sa suffix array of s
 * @param n length of s
 * @param useCS flag for using constant stack
 */
void sa2nsv(const unsigned char * s, unsigned int * sa, int n, bool useCS = false);

/** 
 * compute nsv array from text
 * 
 * @param s input text
 * @param SA suffix array of s
 * @param n length of s
 * @param K alphabet size
 * @param useCS 
 */
void text2nsv(const unsigned char *s, unsigned int *SA,
              unsigned int n, unsigned int K,
              bool  useCS = false);
void text2LMSsa(const unsigned char *s, unsigned int *SA,
                unsigned int n, unsigned int K,
                unsigned int & n1);

void nsv2pnsv(unsigned int n, unsigned int * in_nsv,
              unsigned int * psv, unsigned int * nsv);
