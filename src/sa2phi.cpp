////////////////////////////////////////////////////////////////////////////////
// sa2phi.cpp
////////////////////////////////////////////////////////////////////////////////
#include <vector>
#include <cassert>
#include <iostream>

#include "sa2phi.hpp"
#include "saca-k.hpp"

#define STACK_BITS 16
#define STACK_SIZE (1 << STACK_BITS)
#define STACK_HALF (1 << (STACK_BITS - 1))
#define STACK_MASK ((STACK_SIZE) - 1)


void show(int *A, int n){
  for(int i = 0; i < n; i++){
    std::cout << std::max(-1, A[i]) << ", ";
  }
  std::cout << std::endl;
}


// void sa2LMSsa(const unsigned char * s, int * sa, int n, int & n1){
void sa2LMSsa(const unsigned char * s, unsigned int * sa, unsigned int n, unsigned int & n1){
  n1 = 0;
  // const int n = s.size();
  std::vector<int> Sbkts(CHAR_SIZE, 0); // Sbkts[255] and Sbkte are never used
  std::vector<int> Sbkte(CHAR_SIZE, 0);
  int i;

  // sampling LMS-suffix from SA
  bool cur_type = LTYPE;
  Sbkte[s[n-1]]++;
  for(i = n-2; i >= 0; i--){
    // bkte[s[i]]++;
    if ( s[i] <  s[i+1]) cur_type = STYPE;
    if ( s[i] >  s[i+1]) cur_type = LTYPE;
    // if((s[i] < s[i+1] && (unsigned char) s[i] > (unsigned char) s[i+1]) ||
    //    (s[i] > s[i+1] && (unsigned char) s[i] < (unsigned char) s[i+1])){
    //   std::cout << "[" << (int)s[i] << ", " << (int)s[i+1] << "]=("
    //             << (unsigned char) s[i] << ", " << (unsigned char) s[i+1] << ")"
    //             << std::endl;
    // }
    // if(s[i] < s[i+1]) assert((unsigned char) s[i] < (unsigned char) s[i+1]);
    // if(s[i] > s[i+1]) assert((unsigned char) s[i] > (unsigned char) s[i+1]);
    // assert(0 <= i && i < n);
    if (cur_type == STYPE){ Sbkts[s[i]]++;}
    else {Sbkte[s[i]]++;}
    // sum1++;
    // types[i] = cur_type;
  }
  // std::cout << "sum1=" << sum1
  //           << " sum2=" << sum2 << std::endl;

  unsigned int sum = 0;
  for (i = 0; i < CHAR_SIZE; i++){
    sum += Sbkts[i] + Sbkte[i];
    Sbkts[i] = sum - Sbkts[i]+1;
    Sbkte[i] = sum;
    // Sbkts[i] = sum - Sbkts[i];
    // Sbkte[i] = sum-1;
    // the interval between Sbkts[i] and Sbkte[i] equals to
    // the interval of the suffixes prefixed by i and Stype in SA
    // Note Sbkte[i] may be negative integer
  }
  // std::cout << "n=" << n << " sum=" << sum << std::endl;
  assert(n == sum);
  for(i = 0; i < CHAR_SIZE; i++){
    for (int j = Sbkts[i]; j <= Sbkte[i]; j++){
      // std::cout << "i=" << i << " j=" << j << " Sbkts[i]=" << Sbkts[i]
      //           << " Sbkte[i]=" << Sbkte[i] << std::endl;
      assert(j > 0);
      if (sa[j-1] > 0 && Sbkts[s[sa[j-1]-1]] > j-1){
        // s[sa[j]-1] is Ltype, so sa[j] is LMS-suffix
        sa[n1++] = sa[j-1];
      }
    }
  }
  assert(2*n1 < n);
  for(unsigned int i = n1; i < n; i++){
    sa[i] = EMPTY;
  }
}

// void sa2LMSphi(const unsigned char * s, int * sa, int n, int n1){
void sa2LMSphi(const unsigned char * s, unsigned int * sa, unsigned int n, unsigned int n1){
  // int n = s.size();
  int sa0 = sa[0];
  int i;
  // Now, sa[0..n1-1] must be the suffix array of LMS-suffix

  for(i = n1-1; i > 0; i--){
    sa[2*i] = sa[i];
    sa[2*i-1] = EMPTY;
  }
  int next = EMPTY;
  for(i = n1-1; i >= 0; i--){
    const unsigned int cur = sa[2*i];
    sa[2*i] = EMPTY;
    // sa[2*i] = EMPTY;
    // int next = (i+1 == n1) ? EMPTY : sa[2*(i+1)]; // 2*n1 < nなら条件分岐いらないけど・・・
    // int next = sa[2*(i+1)]; // 2*n1 < nなら条件分岐いらないけど・・・
    // だめだsa[2*n1] = EMPTYとはかぎらない
    // 条件分岐は削除できるけど少し面倒
    assert(cur != EMPTY && cur > 0);
    // const int next = sa[2*(i+1)]; // cur +(cur%2)-1 == nextもありえる
    // sa[2*(i+1)] = EMPTY;

    if (sa[cur] == EMPTY) sa[cur] = next;
    else sa[cur-1] = next;
    // sa[cur + (cur % 2) -1] = next;
    next = cur;

    // sa[2*(i+1)] = EMPTY;
    // if (sa[cur] == EMPTY){ // 実は遇奇の判定で行ける
    //   sa[cur] = next;
    // }else{ // sa[cur] is already set by SA of LMS-suffix
    //   assert(sa[cur-1] == EMPTY);
    //   sa[cur-1] = next;
    // }
  }
  // show(sa, n);
  // sa[0] = EMPTY;

  // シーケンシャルにアクセスし，偶数番目を消去，奇数番目をシフトしたほうがキャッシュ効率的にも良いのかも
  // いや単純にはだめだ．少し考えないと．
  // sa[2*i-1]に値が入っていればその値は本当はsa[2*i]に入れなければいけない．
  // sa[2*i] = EMPTYとする所をsa[2*i-1]に値が入っているのかを示すフラグを入れておけば
  // 次の正しい値を入れる作業がsequentialに行えて速いのかもしれない．
  // やっぱむずけえや，わかんね．

  // rearrange the value which is set at a temporary position.
  // sa[0..n-1]全てを走査するのと
  // sa1[0..n1-1]を走査するのどちらが速いだろう？
  // XXXXXXXXXX できればシーケンシャルにやりたい！！ XXXXXXXXXXXX
  unsigned int cur = sa0;
  do{
    // if (cur % 2 == 0){ // sa[2*i] = EMPTYに設定していなくても大丈夫
    if (sa[cur] == EMPTY){ // cur must be even
      // assert(sa[cur-1] != EMPTY);
      // cur == sa[n1-1] sa[cur] and sa[cur-1] == EMPTY
      assert(cur > 0);
      sa[cur] = sa[cur-1];
      sa[cur-1] = EMPTY;
    }
    // std::cout << cur << ", " << sa[cur] <<  std::endl;
    cur = sa[cur];
  }while (cur != EMPTY);



  // ---ボツ
  // for(i = n1-1; i >= 0; i--){
  //   const int cur = sa[2*i];
  //   assert(cur != EMPTY);
  //   const int next = sa[2*(i+1)]; // cur +(cur%2)-1 == nextもありえる
  //   sa[2*(i+1)] = EMPTY;
  //   sa[cur + (cur % 2) -1] = next;
  //   // --------- test
  //   if (cur % 2 == 0){
  //     if (sa[cur] == EMPTY) sa[cur] = next;
  //     else sa[cur-1] = next;
  //   }else{
  //     if (sa[cur+1] == EMPTY) sa[cur+1] = next;
  //     else sa[cur] = next;
  //   }
  // }
}

void inducePhiL(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                std::vector<unsigned int> & Lbkts, std::vector<unsigned int> & Lbkte){
  // int n = s.size();
  Lbkts[s[n-1]] = n-1; // we assume that s[n] is a sentinel that is smallest character in s.
  Lbkte[s[n-1]] = n-1;
  int phi_idx = sa0;
  // std::cout << EMPTY << std::endl;
  for (int i = 0; i < CHAR_SIZE; i++){
    for (int j = 0; j < 2; j++){
      unsigned int cur;
      if (j == 0) cur = Lbkts[i];
      else cur = phi_idx;
      if (cur == EMPTY) continue;

      int prev = EMPTY;
      // while (cur != EMPTY){
      while (cur != EMPTY && s[cur] == i){ // We have to check s[cur] is i if we don't use S type bucket.
        if (cur > 0 &&  s[cur-1] >=s[cur]){
          const unsigned char c = s[cur-1];
          // std::cout << "i=" << i << " j=" << j
          //           << " lbkts=" << Lbkts[c] << " lbkte=" << Lbkte[c]
          //           << " cur=" << cur << " s[cur-1]=" << (int)c
          //           << " phi[lbkte[c]]="
          //           << ((Lbkte[c] == EMPTY) ? -1 : phi[Lbkte[c]])
          //           << std::endl;
          if (Lbkts[c] == EMPTY) Lbkts[c] = cur-1;
          else{
            assert(Lbkte[c] != EMPTY && phi[Lbkte[c]] == EMPTY);
            phi[Lbkte[c]] = cur-1;
          }
          Lbkte[c] = cur-1;
        }
        unsigned int next = phi[cur];
        assert(cur != next);
        if (j == 0) phi[cur] = prev;
        else  phi[cur] = EMPTY; // s[cur] is LMS-suffix, XXXX[it is really needed?]XXXX
        prev = cur;
        cur = next;
      }
      if (j == 1) {
        // assert(prev == EMPTY || s[prev] == i && s[cur] != i);
        phi_idx = cur;
      }
    }
  }
}


void inducePhiS(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte){
// Sbkts[c] indicate the position in T, where lexicographically smallest S-type suffix of prefix c starts.
// Sbkte[c] indicate the position in T, where lexicographically largest S-type suffix of prefix c starts.
  std::vector<unsigned int> Sbkts(CHAR_SIZE, EMPTY); // Sbkts[255] and Sbkte are never used
  std::vector<unsigned int> Sbkte(CHAR_SIZE, EMPTY);
  assert(phi[sa0] == EMPTY);
  unsigned int prev = EMPTY;
  // int prev = sa0;
  for (int i = CHAR_SIZE-1; i >= 0; i--){
    for (int j = 0; j < 2; j++){
      unsigned int cur;
      if (j == 0) cur = Sbkte[i];
      else cur = Lbkte[i];
      if (cur == EMPTY) continue;

      if (prev != EMPTY) phi[prev] = cur;
      while (cur != EMPTY){
        assert(cur < n);
        // if (cur > 0 && ((j == 0 && (unsigned char) s[cur-1] <= (unsigned char) s[cur]) ||
        //                 ((j == 1) && (unsigned char) s[cur-1] < (unsigned char) s[cur]))){
        if (cur > 0 && (s[cur-1]+j <= s[cur])){
          // if (cur > 0 && type[cur-1] == STYPE){
          // s[cur-1] is S type
          const unsigned char c = s[cur-1];
          if (Sbkte[c] == EMPTY) Sbkte[c] = cur-1;
          else{
            const int prev_S = Sbkts[c];
            //   if(!(prev_S == sa0 || phi[prev_S] == EMPTY)){
            // std::cout << "EMPTY=" << EMPTY << std::endl;
            // std::cout << "i=" << i << " j=" << j
            //           << " cur=" << cur << " n=" << s.size()
            //           << "phi[" << prev_S << "]=" << phi[prev_S]
            //           << std::endl;

            //   }
            assert(prev_S == sa0 || phi[prev_S] == EMPTY);
            phi[prev_S] = cur-1;
          }
          Sbkts[c] = cur-1;
        }
        const int next = phi[cur];
        prev = cur;
        cur = next;
      }
    }
  }
}


void sa2phi(const unsigned char * s, unsigned int * sa, unsigned int n){

  std::vector<unsigned int> Lbkts(CHAR_SIZE, EMPTY);
  std::vector<unsigned int> Lbkte(CHAR_SIZE, EMPTY);
  unsigned int n1;
  int sa0 = sa[0];
  sa2LMSsa(s, sa, n, n1);
  int lmssa0 = sa[0];

  sa2LMSphi(s, sa, n, n1);
  // sa2LMSphi_old(s, sa, sa0);
  std::cout << "finish sa to LMS phi" << std::endl;

  inducePhiL(s, sa, n, lmssa0, Lbkts, Lbkte);
  std::cout << "finish induce Phi L" << std::endl;

  inducePhiS(s, sa, n, lmssa0, Lbkts, Lbkte);
  sa[sa0] = EMPTY; // how define sa[sa0]? EMPTY or last_index?
  std::cout << "finish induce Phi S" << std::endl;
}


// Confirm the definition of sa[sa0] and
// the value of psv_text when it is none value.
// I supppose that sa[sa0] is EMPTY in a former situation,
// and the value is -1 in a latter situation

void nsvTOFromPhi(unsigned int *phi, unsigned int last_idx){
  unsigned int prev = EMPTY;
  unsigned int cur;
  int next;
  cur = last_idx;
  while(cur != EMPTY){
    while(prev != EMPTY && cur < prev){
      // std::cout << cur << ", " << prev <<std::endl;
      assert(prev != EMPTY);
      prev = phi[prev];
    }
    next = phi[cur];
    phi[cur] = prev;
    prev = cur;
    cur = next;
  }
}


void peakElim_ow(int j, int i, unsigned int * nsv, int bot){
  if(j < i){
    peakElim_ow(j, nsv[i], nsv, bot);
  } else {
    const int pre = nsv[j];
    nsv[j] = i;
    if(pre != bot){
      peakElim_ow(pre, j, nsv, bot);
    }
  }
}

void phiL2nsv_CS(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
                 const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte){
  std::vector<unsigned int> Sbkts(CHAR_SIZE, EMPTY); // Sbkts[255] and Sbkte are never used
  std::vector<unsigned int> Sbkte(CHAR_SIZE, EMPTY);
  int * stack = new int[STACK_SIZE + 5];
  // int * stack = new int[STACK_SIZE + 5];
  unsigned int top = 0;
  stack[top] = 0;
  // stack[top] = EMPTY;
  assert(phi[sa0] == EMPTY);
  unsigned int prev = EMPTY;
  for (int i = CHAR_SIZE-1; i >= 0; i--){
    for (int j = 0; j < 2; j++){
      unsigned int cur;
      if (j == 0) cur = Sbkte[i];
      else cur = Lbkte[i];
      if (cur == EMPTY) continue;
      // if (prev != EMPTY) phi[prev] = cur;
      while (cur != EMPTY){
        assert(cur < n);
        // induced sort
        if (cur > 0 && (s[cur-1]+j <= s[cur])){
          const unsigned char c = s[cur-1];
          if (Sbkte[c] == EMPTY) Sbkte[c] = cur-1;
          else{
            const int prev_S = Sbkts[c];
            assert(prev_S == sa0 || phi[prev_S] == EMPTY);
            phi[prev_S] = cur-1;
          }
          Sbkts[c] = cur-1;
        }
        // calculating nsv

        // NSVの計算
        // ここの実装をどうするか？
        // 再帰，while, constant stack?
        while(stack[top] > (cur+1)) top--;
        if ((top & STACK_MASK) == 0){
          if(stack[top] < 0){
            top -= stack[top];
            while (top > (cur+1)) top = phi[top];
            stack[0] = -phi[top];
            stack[1] = top;
            top = 1;
          } else if (top == STACK_SIZE) {
            // Stack is full -- discard half.
            for (int j = STACK_HALF; j <= STACK_SIZE; ++j)
              stack[j - STACK_HALF] = stack[j];
            stack[0] = -stack[0];
            top = STACK_HALF;
          }
        }
        const int next = phi[cur];
        phi[cur] = stack[top] > 0 ? stack[top]-1 : EMPTY;
        top++;
        stack[top] = cur+1;
        prev = cur;
        cur = next;
      }
    }
  }
  delete [] stack;
}
void phiL2nsv(const unsigned char * s, unsigned int * phi, unsigned int n, int sa0,
              const std::vector<unsigned int> & Lbkts, const std::vector<unsigned int> & Lbkte){
// Sbkts[c] indicate the position in T, where lexicographically smallest S-type suffix of prefix c starts.
// Sbkte[c] indicate the position in T, where lexicographically largest S-type suffix of prefix c starts.
  std::vector<unsigned int> Sbkts(CHAR_SIZE, EMPTY); // Sbkts[255] and Sbkte are never used
  std::vector<unsigned int> Sbkte(CHAR_SIZE, EMPTY);
  assert(phi[sa0] == EMPTY);
  unsigned int prev = EMPTY;
  for (int i = CHAR_SIZE-1; i >= 0; i--){
    for (int j = 0; j < 2; j++){
      unsigned int cur;
      if (j == 0) cur = Sbkte[i];
      else cur = Lbkte[i];
      if (cur == EMPTY) continue;
      // if (prev != EMPTY) phi[prev] = cur;
      while (cur != EMPTY){
        assert(cur < n);
        if (cur > 0 && (s[cur-1]+j <= s[cur])){
          const unsigned char c = s[cur-1];
          if (Sbkte[c] == EMPTY) Sbkte[c] = cur-1;
          else{
            const int prev_S = Sbkts[c];
            assert(prev_S == sa0 || phi[prev_S] == EMPTY);
            phi[prev_S] = cur-1;
          }
          Sbkts[c] = cur-1;
        }
        // calculating nsv
        const int next = phi[cur];

        // ここの実装をどうするか？
        // 再帰，while, constant stack?
        phi[cur] = peakElimRight(cur, prev, phi);
        prev = cur;
        cur = next;
      }
    }
  }
}

int peakElimRight(unsigned int j, unsigned int i, const unsigned int * nsv){
  if (i == EMPTY || j > i) return i;
  else return peakElimRight(j, nsv[i], nsv);
}

void sa2nsv(const unsigned char * s, unsigned int * sa, int n, bool useCS){
  std::vector<unsigned int> Lbkts(CHAR_SIZE, EMPTY);
  std::vector<unsigned int> Lbkte(CHAR_SIZE, EMPTY);
  unsigned int n1;
  sa2LMSsa(s, sa, n, n1);
  int lmssa0 = sa[0];
  sa2LMSphi(s, sa, n, n1);
  inducePhiL(s, sa, n, lmssa0, Lbkts, Lbkte);
  if (useCS) phiL2nsv_CS(s, sa, n, lmssa0, Lbkts, Lbkte);
  else phiL2nsv(s, sa, n, lmssa0, Lbkts, Lbkte);
  // inducePhiS(s, sa, lmssa0, Lbkts, Lbkte);
}


void text2nsv(const unsigned char *s, unsigned int *SA,
              unsigned int n, unsigned int K,
              bool useCS) {
  unsigned int n1;
  text2LMSsa(s, SA, n, K, n1);
  std::vector<unsigned int> Lbkts(CHAR_SIZE, EMPTY);
  std::vector<unsigned int> Lbkte(CHAR_SIZE, EMPTY);
  int lmssa0 = SA[1];
  // sa2LMSphi(s, SA, n1);
  // castに注意
  for(unsigned int i = n1; i < n; i++) SA[i] = EMPTY;
  sa2LMSphi(s, SA+1, n-1, n1-1);
  inducePhiL(s, SA+1, n-1, lmssa0, Lbkts, Lbkte);
  if (useCS) phiL2nsv_CS(s, SA+1, n-1, lmssa0, Lbkts, Lbkte);
  else phiL2nsv(s, SA+1, n-1, lmssa0, Lbkts, Lbkte);
  // phiL2nsv_CS(s, SA+1, n-1, lmssa0, Lbkts, Lbkte);
}
void text2LMSsa(const unsigned char *s, unsigned int *SA,
                unsigned int n, unsigned int K,
                unsigned int & n1) {
  unsigned int i;
  unsigned int *bkt=NULL;
  unsigned int m = n;
  // stage 1: reduce the problem by at least 1/2.

  bkt=(unsigned int *)malloc(sizeof(int)*K);
  putSubstr0(SA, s, bkt, n, K);
  induceSAl0(SA, s, bkt, n, K, false);
  induceSAs0(SA, s, bkt, n, K, false);

  // now, all the LMS-substrings are sorted and 
  //   stored sparsely in SA.

  // compact all the sorted substrings into
  //   the first n1 items of SA.
  // 2*n1 must be not larger than n.
  // unsigned int n1=0;
  n1 = 0;
  unsigned int level = 0;
  for(i=0; i<n; i++) 
    if((!level&&SA[i]>0) || (level&&((int *)SA)[i]>0))
      SA[n1++]=SA[i];

  unsigned  *SA1=SA, *s1=SA+m-n1;
  unsigned int name_ctr;
  name_ctr=nameSubstr(SA,s,s1,n,m,n1,level);

  // stage 2: solve the reduced problem.

  // recurse if names are not yet unique.
  if(name_ctr<n1)
    SACA_K((unsigned char *)s1, SA1, 
           n1, 0, m-n1, level+1);
  else // get the suffix array of s1 directly.
    for(i=0; i<n1; i++) SA1[s1[i]]=i;

  // stage 3: induce SA(S) from SA(S1).

  getSAlms(SA, s, s1, n, n1, level);
  // Now, SA[1..n1-1] is suffix array of LMS-suffixes
  free(bkt);
}

void nsv2pnsv(unsigned int n, unsigned int * in_nsv,
              unsigned int * psv, unsigned int * nsv){
  unsigned int i;
  unsigned int imax = 0;
  nsv[0] = EMPTY;
  psv[0] = EMPTY;
  in_nsv[0] = EMPTY;
  for(i = 1; i < n; i++){
    // if (in_nsv[i]==EMPTY) std::cout << i << "empty";
    nsv[i] = in_nsv[i];
    if (nsv[i] != EMPTY){
      psv[i] = in_nsv[nsv[i]];
      in_nsv[nsv[i]] = i;
    }else{
      psv[i] = imax;
      imax = i;
    }
    in_nsv[i] = psv[i];
  }
  // now in_nsv is phi array
}
