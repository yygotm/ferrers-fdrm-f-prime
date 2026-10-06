// Simulated annealing for an [F',6,7]_q FDRM code over GF(q), q in {2,3,4,5,7}, F' = [1,2,3,4,8,8,8,8,8].
// Same parametrisation as sa2.c: basis k has a 1 at the k-th dot of columns 0..2, the other 44 dots are free.
// Codewords are enumerated projectively (first nonzero coefficient = 1). Cost = sum max(0, 7-rank).
// GF(4) encoding: 0,1,2=w,3=w^2=w+1.  usage: saq q seed steps T0
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static const int C[9] = {1,2,3,4,8,8,8,8,8};
#define R 8
#define NC 9
#define K 6
static int q;
static uint8_t ADD[8][8], MUL[8][8], INV[8], NEG[8];
static uint8_t B[K][R][NC];
static int freeR[64], freeC[64], nfree;
static int ncw; static uint8_t (*coef)[K]; static uint8_t (*CW)[R][NC]; static int *rk;

static uint64_t rng = 88172645463325252ULL;
static inline uint64_t xr(void){ rng ^= rng<<13; rng ^= rng>>7; rng ^= rng<<17; return rng; }

static void setup_field(void){
  for(int a=0;a<q;a++) for(int b=0;b<q;b++){
    if(q==4){ int r=0; for(int i=0;i<2;i++) if(b>>i&1) r^=a<<i; if(r&4) r^=7; ADD[a][b]=a^b; MUL[a][b]=r; }
    else { ADD[a][b]=(a+b)%q; MUL[a][b]=(a*b)%q; }
  }
  for(int a=0;a<q;a++) for(int b=0;b<q;b++){ if(MUL[a][b]==1) INV[a]=b; if(ADD[a][b]==0) NEG[a]=b; }
}

static int rankq(uint8_t m[R][NC]){
  uint8_t a[R][NC]; memcpy(a,m,sizeof a);
  int r=0;
  for(int c=0;c<NC && r<R;c++){
    int p=-1; for(int i=r;i<R;i++) if(a[i][c]){p=i;break;}
    if(p<0) continue;
    if(p!=r){ uint8_t t[NC]; memcpy(t,a[p],NC); memcpy(a[p],a[r],NC); memcpy(a[r],t,NC); }
    uint8_t s=INV[a[r][c]]; for(int j=c;j<NC;j++) a[r][j]=MUL[s][a[r][j]];
    for(int i=r+1;i<R;i++) if(a[i][c]){ uint8_t f=NEG[a[i][c]]; for(int j=c;j<NC;j++) a[i][j]=ADD[a[i][j]][MUL[f][a[r][j]]]; }
    r++;
  }
  return r;
}

static int pen(int r){ return r<7? 7-r:0; }

int main(int argc,char**argv){
  q=atoi(argv[1]); uint64_t seed=strtoull(argv[2],0,10); long steps=atol(argv[3]); double T0=atof(argv[4]);
  setup_field();
  rng ^= seed*0x9E3779B97F4A7C15ULL; for(int i=0;i<10;i++) xr();
  int fr[6]={0,0,1,0,1,2}, fc[6]={0,1,1,2,2,2};
  memset(B,0,sizeof B);
  for(int k=0;k<K;k++) B[k][fr[k]][fc[k]]=1;
  for(int c=3;c<NC;c++) for(int r=0;r<C[c];r++){ freeR[nfree]=r; freeC[nfree]=c; nfree++; }
  for(int k=0;k<K;k++) for(int f=0;f<nfree;f++) B[k][freeR[f]][freeC[f]]=xr()%q;
  // projective coefficient vectors
  int tot=1; for(int k=0;k<K;k++) tot*=q;
  coef=malloc(tot*K); ncw=0;
  for(int x=1;x<tot;x++){ uint8_t v[K]; int y=x; for(int k=0;k<K;k++){ v[k]=y%q; y/=q; }
    int fnz=0; while(!v[fnz]) fnz++; if(v[fnz]!=1) continue; memcpy(coef[ncw++],v,K); }
  CW=malloc(sizeof(*CW)*ncw); rk=malloc(sizeof(int)*ncw);
  int cur=0;
  for(int w=0;w<ncw;w++){ memset(CW[w],0,sizeof CW[w]);
    for(int k=0;k<K;k++) if(coef[w][k]) for(int i=0;i<R;i++) for(int j=0;j<NC;j++) CW[w][i][j]=ADD[CW[w][i][j]][MUL[coef[w][k]][B[k][i][j]]];
    rk[w]=rankq(CW[w]); cur+=pen(rk[w]); }
  int best=cur;
  int *newrk=malloc(sizeof(int)*ncw);
  for(long s=0;s<steps;s++){
    double T=T0*(1.0-(double)s/steps)+1e-3;
    int k=xr()%K, f=xr()%nfree, i=freeR[f], j=freeC[f];
    uint8_t old=B[k][i][j], nv=(old+1+xr()%(q-1))%q;   // new value != old
    uint8_t delta=ADD[nv][NEG[old]];
    int nc=cur;
    for(int w=0;w<ncw;w++) if(coef[w][k]){
      uint8_t sv=CW[w][i][j]; CW[w][i][j]=ADD[sv][MUL[coef[w][k]][delta]];
      newrk[w]=rankq(CW[w]); nc+=pen(newrk[w])-pen(rk[w]); CW[w][i][j]=sv; }
    if(nc<=cur || exp((cur-nc)/T) > (xr()%1000000)/1e6){
      B[k][i][j]=nv; cur=nc;
      for(int w=0;w<ncw;w++) if(coef[w][k]){ CW[w][i][j]=ADD[CW[w][i][j]][MUL[coef[w][k]][delta]]; rk[w]=newrk[w]; }
      if(cur<best){ best=cur; fprintf(stderr,"q=%d step %ld cost %d\n",q,s,cur); }
    }
    if(cur==0){
      printf("FOUND q=%d seed=%llu step=%ld\n",q,(unsigned long long)seed,s);
      for(int kk=0;kk<K;kk++){ for(int r=0;r<R;r++){ for(int c=0;c<NC;c++) putchar('0'+B[kk][r][c]); putchar(r<R-1?' ':'\n'); } }
      return 0;
    }
  }
  printf("NOTFOUND q=%d seed=%llu best=%d\n",q,(unsigned long long)seed,best);
  return 1;
}
