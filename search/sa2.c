// Simulated annealing for a binary [F',6,7] FDRM code, F' = [1,2,3,4,8,8,8,8,8] (8 rows x 9 cols, top-aligned).
// The code is the graph of a linear map from the 6 dots of columns 0..2 (forced injective by the
// Etzion-Silberstein bound at t=0) to the remaining 44 dots. Cost = sum over nonzero codewords of max(0, 7-rank).
// usage: sa2 seed steps T0
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

static const int C[9] = {1,2,3,4,8,8,8,8,8};
#define R 8
#define K 6
static uint16_t B[K][R];          // basis matrices, row bitmasks over 9 columns
static int freeR[64], freeC[64], nfree;

static uint64_t rng = 88172645463325252ULL;
static inline uint64_t xr(void){ rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17; return rng; }

static int rank8(const uint16_t *m){
  uint16_t a[R]; for(int i=0;i<R;i++) a[i]=m[i];
  int r=0;
  for(int c=0;c<9 && r<R;c++){
    int p=-1; for(int i=r;i<R;i++) if(a[i]>>c&1){p=i;break;}
    if(p<0) continue;
    uint16_t t=a[p]; a[p]=a[r]; a[r]=t;
    for(int i=0;i<R;i++) if(i!=r && (a[i]>>c&1)) a[i]^=a[r];
    r++;
  }
  return r;
}

static int cost(int *minrank){
  int s=0, mn=99; uint16_t m[R];
  for(int x=1;x<(1<<K);x++){
    for(int i=0;i<R;i++){ uint16_t v=0; for(int k=0;k<K;k++) if(x>>k&1) v^=B[k][i]; m[i]=v; }
    int rk=rank8(m); if(rk<mn) mn=rk; if(rk<7) s+=7-rk;
  }
  if(minrank) *minrank=mn;
  return s;
}

int main(int argc,char**argv){
  uint64_t seed = argc>1? strtoull(argv[1],0,10):1; long steps = argc>2? atol(argv[2]):10000000; double T0 = argc>3? atof(argv[3]):1.0;
  rng ^= seed*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) xr();
  // fixed part: dots of columns 0..2 in order (r,c): (0,0),(0,1),(1,1),(0,2),(1,2),(2,2) -> basis k has a 1 at the k-th
  int fr[6]={0,0,1,0,1,2}, fc[6]={0,1,1,2,2,2};
  for(int k=0;k<K;k++){ for(int i=0;i<R;i++) B[k][i]=0; B[k][fr[k]] |= 1u<<fc[k]; }
  for(int c=3;c<9;c++) for(int r=0;r<C[c];r++){ freeR[nfree]=r; freeC[nfree]=c; nfree++; }
  for(int k=0;k<K;k++) for(int f=0;f<nfree;f++) if(xr()&1) B[k][freeR[f]] ^= 1u<<freeC[f];
  int mr, cur=cost(&mr), best=cur;
  for(long s=0;s<steps;s++){
    double T = T0*(1.0 - (double)s/steps) + 1e-3;
    int k=xr()%K, f=xr()%nfree;
    B[k][freeR[f]] ^= 1u<<freeC[f];
    int nc=cost(&mr);
    if(nc<=cur || exp((cur-nc)/T) > (xr()%1000000)/1e6){ cur=nc; if(cur<best){best=cur; fprintf(stderr,"step %ld cost %d\n",s,cur);} }
    else B[k][freeR[f]] ^= 1u<<freeC[f];
    if(cur==0){
      cost(&mr);
      printf("FOUND seed=%llu step=%ld minrank=%d\n",(unsigned long long)seed,s,mr);
      for(int kk=0;kk<K;kk++){ for(int i=0;i<R;i++){ for(int c=0;c<9;c++) putchar(B[kk][i]>>c&1?'1':'0'); putchar(i<R-1?' ':'\n'); } }
      return 0;
    }
  }
  printf("NOTFOUND seed=%llu best=%d\n",(unsigned long long)seed,best);
  return 1;
}
