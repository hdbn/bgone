#pragma once
// Author: Ge Nong,  Email: issng@mail.sysu.edu.cn
// Department of Computer Science, Sun Yat-sen University, 
// Guangzhou, China
// Date: December 24, 2012
//
// This is the demo source code for the algorithm SACA-K presented in this article:
// G. Nong, Practical Linear-Time O(1)-Workspace Suffix Sorting for Constant Alphabets, 
// ACM Transactions on Information Systems, Scheduled to Appear in July 2013.
// A draft for this article can be retrieved from http://code.google.com/p/ge-nong/.
#include <stdlib.h>

// // set only the highest bit as 1, i.e. 1000...
// const unsigned int EMP=((unsigned int)1)<<(sizeof(unsigned int)*8-1); 

// // get s[i] at a certain level
// #define chr(i) ((level==0)?((unsigned char *)s)[i]:((int *)s)[i])

void getBuckets(const unsigned char *s, 
                unsigned int *bkt, unsigned int n,
                unsigned int K, bool end);

void putSuffix0(unsigned int *SA, 
                const unsigned char *s, unsigned int *bkt, 
                unsigned int n, unsigned int K, int n1) ;

void induceSAl0(unsigned int *SA,
                const unsigned char *s, unsigned int *bkt,
                unsigned int n, unsigned int K, bool suffix) ;
void induceSAs0(unsigned int *SA,
                const unsigned char *s, unsigned int *bkt,
                unsigned int n, unsigned int K, bool suffix) ;

void putSubstr0(unsigned int *SA,
                const unsigned char *s, unsigned int *bkt,
                unsigned int n, unsigned int K) ;
void putSuffix1(int *SA, int *s, int n1) ;
void induceSAl1(int *SA, int *s, 
                int n, bool suffix) ;

void induceSAs1(int *SA, int *s, 
                int n, bool suffix) ;
void putSubstr1(int *SA, int *s, int n);
unsigned int getLengthOfLMS(const unsigned char *s, 
                            unsigned int n, int level, unsigned int x);

unsigned int nameSubstr(unsigned int *SA, 
                        const unsigned char *s, unsigned int *s1, unsigned int n, 
                        unsigned int m, unsigned int n1, int level) ;

void getSAlms(unsigned int *SA, 
              const unsigned char *s, 
              unsigned int *s1, unsigned int n, 
              unsigned int n1, int level );


void SACA_K(const unsigned char *s, unsigned int *SA,
            unsigned int n, unsigned int K,
            unsigned int m, int level) ;
