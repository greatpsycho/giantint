/**************************************************************

    ellmul.c

    High-performance elliptic multiplication
    modulo Fermat numbers.

    Change history:

    28 Jul 99   ARP - Added conditional variable declarations
    28 Jul 99   REC - Made JSP's separate-file ellmul0.c
    21 Jul 99   JSP - Creation and initial debug/optimize.

    c. 1999 Perfectly Scientific, Inc.
    All Rights Reserved.

**************************************************************/
 
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "giants.h"
 
#define s_modg(N, g) fer_mod(Q, g, N)
 
typedef struct {
   double R, I;
} cplx;
 
#define BASE 65536
#define OUR_PI ((double) 3.141592653589793238462643383279502884197)
 
/*----------------- global variables ------------------------*/
 
cplx *fftcache, *x, *roots, *dwt, *temp,
     *cacheptr, *temptr;
cplx carries[4];
void (*bottom_out)(cplx *c);
double magic[4];
 
extern giant xs, zs, N;                     /* from fermat.c */
extern int Q;
 
/*--------------- prototypes from fermat.c ------------------*/
 
void        ell_even( giant x1, giant z1, giant x2, giant z2,
                      giant An, giant Ad, giant N );
void        ell_odd( giant x1, giant z1, giant x2, giant z2,
                     giant xor, giant zor, giant N );
 
/*--------------- local prototypes --------------------------*/
 
void init_ellmul( long size );
void ell_bit_set( long size );
void ell_bit_clear( long size );
void giant2packed( giant a, cplx *b, long j, long size );
void packed2giant( giant a, cplx *b, long j, long size );
void single_fft2( cplx *c0, cplx *c1 );
void vector_fft2( cplx *c0, cplx *c1 );
void single_splitrad_DIF( cplx *c, cplx *w, long size );
void vector_splitrad_DIF( cplx *c, cplx *w, long size );
void single_splitrad_DIT( cplx *c, cplx *w, long size );
void vector_splitrad_DIT( cplx *c, cplx *w, long size );
void base1(cplx *c);
extern void base2(cplx *c);
extern void base3(cplx *c);
void twist( cplx *c, long size, long vec );
void untwist( cplx *c, long size, long vec );
void vector_fft_ifft( cplx *c, cplx *w, long size );
void single_fft( cplx *c, cplx *w, long size );
void vector_fft( cplx *c, cplx *w, long size );
void single_ifft( cplx *c, cplx *w, long size );
void release_carries( cplx *c, long size, double scale, long vec );
void last_carry( cplx *c, long vec );
 
/*-----------------------------------------------------------*/

void ell_mul( giant xx, giant zz, int n,
              giant An, giant Ad, giant  N) {
 
   unsigned long c = (unsigned long)(0x80000000), i, size = Q/16;
 
   if (n==1)
      return;
   if (n==2) {
      ell_even(xx, zz, xx, zz, An, Ad, N);
      return;
   }
 
   ell_even(xx, zz, xs, zs, An, Ad, N);
   while ((c&n) == 0) c >>= 1;
   c>>=1;
 
   giant2packed( xs, x, 4, size );    /* interleave the four variables */
   giant2packed( zs, x+1, 4, size );
   giant2packed( xx, x+2, 4, size );
   giant2packed( zz, x+3, 4, size );
 
   giant2packed( An, fftcache, 4, size );   /* and the four constants  */
   giant2packed( Ad, fftcache+1, 4, size ); /* zorg=zz, xorg=xx for now*/
   giant2packed( zz, fftcache+2, 4, size );
   giant2packed( xx, fftcache+3, 4, size );
 
   size = size / 2;                    /* size is now size of cplx FFT */
 
   twist( fftcache, size, 4 );
   vector_fft( fftcache, roots, size );     /* cache FFTs of constants */
 
   for(i=0; i<4*size; i+=4) {                   /* pre-scale the cache */
      fftcache[i+1].R *= 4; fftcache[i+1].I *= 4;
      fftcache[i+2].R *= 4; fftcache[i+2].I *= 4;
      fftcache[i+3].R *= 4; fftcache[i+3].I *= 4;
   }
 
   do {
      if(c&n) {
         ell_bit_set(size);
      }
      else {
         ell_bit_clear(size);
      }
      c >>= 1;
   } while(c);
 
   size = 2*size;
   packed2giant( xx, x+2, 4, size );
   packed2giant( zz, x+3, 4, size );       /* salvage the answer */
}
 
/*---------------------------------------------------------------------*/

void special_mulg( giant a, giant b ) {
 
   long size = Q/16, mul_size = size/2, i;
   cplx *aptr = x+mul_size, *bptr = x;
   double f0,f1,f2,f3,f4,f5,scale = 1.0/mul_size;
 
   s_modg(N, a);
   s_modg(N, b);
 
   giant2packed( a, aptr, 1, size );
   twist( aptr, mul_size, 1 );
   single_fft( aptr, roots, mul_size );
 
   giant2packed( b, bptr, 1, size );
   twist( bptr, mul_size, 1 );
   single_fft( bptr, roots, mul_size );
 
   for(i=0; i<mul_size; i++) {
      f0 = bptr[i].R;   f1 = aptr[i].R;
      f2 = bptr[i].I;   f3 = aptr[i].I;
      f4 = f0 * f1;     f5 = f2 * f3;
      f0 = f0 * f3;     f1 = f1 * f2;
      f4 -= f5;         f0 += f1;
      bptr[i].R = f4;   bptr[i].I = f0;
   }
 
   single_ifft( bptr, roots, mul_size );
   untwist( bptr, mul_size, 1 );
   release_carries( bptr, mul_size, scale, 1 );
   last_carry(bptr, 1);
   packed2giant( b, bptr, 1, size );
}
 
/*---------------------------------------------------------------------*/

void special_squareg( giant b ) {
 
   long size = Q/16, mul_size = size/2, i;
   cplx  *bptr = x;
   double f0,f1,f4, scale = 1.0/mul_size;
 
   s_modg(N, b);
 
   giant2packed( b, bptr, 1, size );
   twist( bptr, mul_size, 1 );
   single_fft( bptr, roots, mul_size );
 
   for(i=0; i<mul_size; i++) {
      f0 = bptr[i].R;   f1 = bptr[i].I;
      f4 = f0 * f1;     f0 = f0 * f0;    f1 = f1 * f1;
      f0 -= f1;         f4 += f4;
      bptr[i].R = f0;   bptr[i].I = f4;
   }
 
   single_ifft( bptr, roots, mul_size );
   untwist( bptr, mul_size, 1 );
   release_carries( bptr, mul_size, scale, 1 );
   last_carry(bptr, 1);
   packed2giant( b, bptr, 1, size );
}
 
/*---------------------------------------------------------------------*/

void ell_bit_set( long size ) {
 
   /* this mess does ell_even() and ell_odd() simultaneously,
      using 25 FFTs (instead of the 36 or so that would be
      needed if each mult or square were done individually) */
 
   long i, j;
   double scale = 1.0/size;
 
   bottom_out = base1;                 
   twist(x,size,4);                    /*  x1        (x1+z1)^2     */
   vector_fft_ifft(x, roots, size);    /*  z1 ->     (x1-z1)^2     */
   untwist(x, size, 4);                /*  x2      x1*x2 - z1*z2   */
   release_carries(x, size, scale, 4); /*  z2      x1*z2 - z1*x2   */
   last_carry(x, 4);
 
   for(i=0,j=0; i<size; i++,j+=4) {                
      temp[i].R = x[j].R - x[j+1].R;          /* temp <- 4*x1*z1 */
      temp[i].I = x[j].I - x[j+1].I;     
   }                                     
   twist(temp,size,1);
   single_fft(temp,roots,size);          /* do forward DWT for temp */
   temptr = temp;
 
   cacheptr = fftcache;
   bottom_out = base2;
   twist(x, size, 4);                  /*  x1         x1*z1           */
   vector_fft_ifft(x, roots, size);    /*  z1 -> 4*Ad*z1 + An*(x1-z1) */
   untwist(x, size, 4);                /*  x2          x2^2           */
   release_carries(x, size, scale, 4); /*  z2          z2^2           */
   last_carry(x, 4);
 
   cacheptr = fftcache;             
   bottom_out = base3;
   twist(x,size,4);                    /*  x1         x1*4*Ad         */
   vector_fft_ifft(x, roots, size);    /*  z1 ->      z1*temp         */
   untwist(x, size, 4);                /*  x2        x2*4*zorg        */
   release_carries(x, size, scale, 4); /*  z2        z2*4*xorg        */
   last_carry(x, 4);
}
 
/*---------------------------------------------------------------------*/

void ell_bit_clear( long size ) {
   
   long i;
   double f0,f1,f2,f3;
 
   for(i=0; i<4*size; i+=4) {
      f0 = x[i].R;    f1 = x[i].I;     /* switch the x's and the z's */
      f2 = x[i+2].R;  f3 = x[i+2].I;
      x[i].R = f2;    x[i].I = f3;
      x[i+2].R = f0;  x[i+2].I = f1;
      f0 = x[i+1].R;  f1 = x[i+1].I;
      f2 = x[i+3].R;  f3 = x[i+3].I;
      x[i+1].R = f2;  x[i+1].I = f3;
      x[i+3].R = f0;  x[i+3].I = f1;
   }
 
   ell_bit_set(size);
 
   for(i=0; i<4*size; i+=4) {
      f0 = x[i].R;    f1 = x[i].I;               /* switch them back */
      f2 = x[i+2].R;  f3 = x[i+2].I;
      x[i].R = f2;    x[i].I = f3;
      x[i+2].R = f0;  x[i+2].I = f1;
      f0 = x[i+1].R;  f1 = x[i+1].I;
      f2 = x[i+3].R;  f3 = x[i+3].I;
      x[i+1].R = f2;  x[i+1].I = f3;
      x[i+3].R = f0;  x[i+3].I = f1;
   }
}
 
/*---------------------------------------------------------------------*/

void init_ellmul( long size ) {

double ieee_round = 6755399441055744.0;
double intel_round = 6755399441055744.0 * 2048.0;
double round_test = 2.7;
double round_correct = 3.0;
 
   /* size is the number of shorts needed to store one
      Fermat residue; the mallocs below assume half that
      number of cplx elements. The factor of 4 in x[]
      and scratch[] is in anticipation of interleaved
      structures */
 
   cplx *temptr;
   double w0, round_me;
   long i, j, inc;
 
   size = size/2;
   w0 = 2*OUR_PI/size;
 
   fftcache = (cplx *)malloc( size * 4 * sizeof(cplx) );
   x = (cplx *)malloc( size * 4 * sizeof(cplx) );
   roots = (cplx *)malloc( size * sizeof(cplx) );
   dwt = (cplx *)malloc( size * sizeof(cplx) );
   temp = (cplx *)malloc( size * sizeof(cplx) );
 
   if(!fftcache || !x || !roots || !dwt || !temp) {
      printf("memory allocation failed\n");
      exit(-1);
   }
 
   for(i=0; i<size; i++) {
      dwt[i].R = cos(.25 * i * w0);
      dwt[i].I = sin(.25 * i * w0);
   }
 
   temptr = roots;
   inc = 1;
   while( size>=4 ) {
      for(j=0; j<size/4; j++) {
         temptr->R = cos(w0 * j * inc);
         temptr->I = -sin(w0 * j * inc);     temptr++;
         temptr->R = cos(w0 * 3 * j * inc);
         temptr->I = -sin(w0 * 3 * j * inc); temptr++;
      }
      size = size / 2;
      inc = inc * 2;
   }
 
   magic[0] = intel_round;          /* figure out dynamically the */
   magic[1] = magic[0];             /* size of an FPU register    */
   round_me = round_test;
   round_me += magic[0];
   round_me -= magic[1];
   if( round_me == round_correct ) {
/*     printf("Using Intel-size 64-bit rounding\n"); */
   }
   else {
      magic[0] = ieee_round;
      magic[1] = magic[0];
      round_me = round_test;
      round_me += magic[0];
      round_me -= magic[1];
      if( round_me != round_correct ) {
/*        printf("No way to round numbers correctly. Bailing out\n"); */
         exit(-1);
      }
/*    printf("Using IEEE 53-bit rounding\n"); */
   }
 
   magic[2] = 1.0/BASE;
   magic[3] = BASE;
 
   for(i=0;i<4;i++) {
      carries[i].R = carries[i].I = 0;
   }
}
 
/*---------------------------------------------------------------------*/

void giant2packed( giant a, cplx *b, long j, long size ) {
 
   /* Transfers the giant "a" into an array of cplx's
      b[], with stride j. Each digit of "a" is converted
      to balanced representation, and stored for use with
      an all-complex DWT. size is the number of shorts in a;
      it's assumed that b has at least half that number of
      cplx elements. */
 
   long borrow = 0, i, temp;
   cplx *btmp = b;
 
   for(i=abs(a->sign); i<size; i++)  a->n[i] = 0;  /* zero unused bits */
 
   for(i=0; i<size/2; i++,b+=j) {
      temp = (long)(a->n[i]) + borrow;
      if( temp >= BASE/2 ) {
         b[0].R = (double)(temp - BASE);
         borrow = 1;
      }
      else {
         b[0].R = (double)temp;
         borrow = 0;
      }
   }
 
   b = btmp;
   for(i=0; i<size/2; i++,b+=j) {
      temp = (long)(a->n[i+size/2]) + borrow;
      if( temp >= BASE/2 ) {
         b[0].I = (double)(temp - BASE);
         borrow = 1;
      }
      else {
         b[0].I = (double)temp;
         borrow = 0;
      }
   }
 
   btmp[0].R -= borrow;
 
   if(a->sign > size)                          /* wrap around possible */
      btmp[0].R -= (double)(a->n[a->sign-1]);    /* previous carry out */
}
 
/*---------------------------------------------------------------------*/

void packed2giant( giant a, cplx *b, long j, long size ) {
 
   /* performs the reverse of giant2packed(), converting
      to standard representation in the process. Unlike that
      code, it's critical that every digit be normalized
      correctly. b[] is assumed to be completely in balanced
      representation, so that all the elements of b[] lie 
      strictly in [-(BASE/2-1), BASE/2-1]  */
 
   long borrow = 0, i, temp;
   cplx *btmp = b;
 
   for(i=0; i<size/2; i++,b+=j) {
      temp = (long)b[0].R + borrow;
      borrow = 0;
      while( temp<0 ) {
         temp += BASE;
         borrow--;
      }
      a->n[i] = (unsigned short)temp;
   }
 
   b = btmp;
   for(i=0; i<size/2; i++,b+=j) {
      temp = (long)b[0].I + borrow;
      borrow = 0;
      while( temp<0 ) {
         temp += BASE;
         borrow--;
      }
      a->n[i+size/2] = (unsigned short)temp;
   }
 
   while( size>=0 && a->n[size-1]==0 )      /* strip off leading zeros */
      size--;                               /* giants.c requires this! */
   a->sign = size;
   iaddg( -borrow, a );                          /* wrap borrow around */
}
 
/*---------------------------------------------------------------------*/

void vector_fft_ifft( cplx *c, cplx *w, long size ) {
 
   /* fully recursive FFT multiply. This does four FFTs
      interleaved together, then performs some sort of
      convolution operation (whatever function *bottom_out
      points to), then the four corresponding IFFTs. This
      is a standard complex split-radix; a more thorough
      implementation would also code radix-4 transforms
      so the recursion bottoms out earlier and saves function
      calls.    */
 
   if( size==2 ) {
      vector_fft2(c, c+4);
 
      (*bottom_out)(c);
      (*bottom_out)(c+4);
 
      vector_fft2(c, c+4);
   }
   else if ( size==4 ) {
      vector_splitrad_DIF( c, w, size );
      vector_fft2(c, c+4);
 
      (*bottom_out)(c);
      (*bottom_out)(c+4);
      (*bottom_out)(c+8);
      (*bottom_out)(c+12);
 
      vector_fft2(c, c+4);
      vector_splitrad_DIT( c, w, size );
   }
   else {
      vector_splitrad_DIF( c, w, size );
      vector_fft_ifft(c, w+size/2, size/2);
      vector_fft_ifft(c+4*size/2, w+3*size/4, size/4);
      vector_fft_ifft(c+4*3*size/4, w+3*size/4, size/4);
      vector_splitrad_DIT( c, w, size );
   }
}
 
/*---------------------------------------------------------------------*/

void single_fft( cplx *c, cplx *w, long size ) {
 
   /* performs a single forward split-radix FFT */
 
   if( size==2 ) {
      single_fft2(c, c+1);
   }
   else if ( size==4 ) {
      single_splitrad_DIF( c, w, size );
      single_fft2(c, c+1);
   }
   else {
      single_splitrad_DIF( c, w, size );
      single_fft(c, w+size/2, size/2);
      single_fft(c+size/2, w+3*size/4, size/4);
      single_fft(c+3*size/4, w+3*size/4, size/4);
   }
}
 
/*---------------------------------------------------------------------*/

void single_ifft( cplx *c, cplx *w, long size ) {
 
   /* performs a single inverse split-radix FFT */
 
   if( size==2 ) {
      single_fft2(c, c+1);
   }
   else if ( size==4 ) {
      single_fft2(c, c+1);
      single_splitrad_DIT( c, w, size );
   }
   else {
      single_ifft(c, w+size/2, size/2);
      single_ifft(c+size/2, w+3*size/4, size/4);
      single_ifft(c+3*size/4, w+3*size/4, size/4);
      single_splitrad_DIT( c, w, size );
   }
}
 
/*---------------------------------------------------------------------*/

void vector_fft( cplx *c, cplx *w, long size ) {
 
   /* four simultaneous forward split-radix FFTs */
 
   if( size==2 ) {
      vector_fft2(c, c+4);
   }
   else if ( size==4 ) {
      vector_splitrad_DIF( c, w, size );
      vector_fft2(c, c+4);
   }
   else {
      vector_splitrad_DIF( c, w, size );
      vector_fft(c, w+size/2, size/2);
      vector_fft(c+4*size/2, w+3*size/4, size/4);
      vector_fft(c+4*3*size/4, w+3*size/4, size/4);
   }
}

/*---------------------------------------------------------------------*/
 
      /*-------------------------------------------------*/
      /* These are the lowest level routines. There's no */
      /* assembly language here, but the C is written so */
      /* that certain compilers produce efficient code.  */
      /* The default is for intel-based gcc, and makes   */
      /* efficient use of the Pentium FPU (the DIT code  */
      /* much more than the DIF code). Everything in the */
      /* RISC section produces efficient code for Sun cc */
      /* on an UltraSPARC, and more aggressive compilers */
      /* will do even better. BE CAREFUL; the two dif-   */
      /* ferent coding styles are much slower when used  */
      /* on the wrong architecture.                      */
      /*                                                 */
      /* Both are also a terrible mess. Abandon all hope */
      /* ye who press enter here.                        */
      /*-------------------------------------------------*/
 
#define RISC_DIF(q)                                                   \
      f10 = w[0].R;   f11 = w[0].I;                                   \
      f12 = w[1].R;   f13 = w[1].I;                                   \
      for(i=0;i<q;i++) {                                              \
         f0 = c0[i].R;   f1 = c0[dist+i].R;  f2=f0+f1; f6=f0-f1;      \
         f0 = c0[i].I;   f1 = c0[dist+i].I;  f3=f0+f1; f7=f0-f1;      \
         f0 = c1[i].I;   f1 = c1[dist+i].I;  f5=f0+f1; f8=f0-f1;      \
         f0 = c1[i].R;   f1 = c1[dist+i].R;  f4=f0+f1; f9=f1-f0;      \
         c0[i].R = f2;   f2=f6+f8;         c0[i].I = f3;   f3=f7+f9;  \
         c1[i].R = f4;   f4=f6-f8;         c1[i].I = f5;   f5=f7-f9;  \
         f6 = f2 * f10;        f7 = f3 * f11;                         \
         f8 = f4 * f12;        f9 = f5 * f13;                         \
         f2 = f2 * f11;        f3 = f3 * f10;    f6 = f6 - f7;        \
         f4 = f4 * f13;        f5 = f5 * f12;    f8 = f8 - f9;        \
         f2 = f2 + f3;        c0[dist+i].R = f6;                      \
         f4 = f4 + f5;        c1[dist+i].R = f8;                      \
         c0[dist+i].I = f2;                                           \
         c1[dist+i].I = f4;                                           \
      }
 
#define PENTIUM_DIF(q)                                                     \
        f0 = c0[q].R;      f1 = c0[dist+q].R;      f0 = f0 - f1;           \
        f1 = c0[q].I;      f2 = c0[dist+q].I;      f1 = f1 - f2;           \
        f2 = c1[q].I;      f3 = c1[dist+q].I;      f2 = f2 - f3;           \
        f3 = c1[q].R;      f4 = c1[dist+q].R;      f3 = f4 - f3;           \
        f4 = f0;           f0 = f0 + f2;         f2 = f4 - f2;             \
        f4 = f1;           f1 = f1 + f3;         f3 = f4 - f3;             \
        f4 = f0;           f5 = w[0].I; f4 = f4 * f5;    f0 = f0 * w[0].R; \
        f5 = f1;           f6 = w[0].R; f5 = f5 * f6;    f1 = f1 * w[0].I; \
        f6 = c0[q].R;      f6 = f6 + c0[dist+q].R; f4 = f4 + f5;           \
        f5 = c0[q].I;      f5 = f5 + c0[dist+q].I; f0 = f0 - f1;           \
        c0[q].R = f6;      f6 = f2;              f6 = f6 * w[1].R;         \
        c0[q].I = f5;      f5 = f3;              f5 = f5 * w[1].I;         \
        c0[dist+q].R = f0;   f2 = f2 * w[1].I;                             \
        c0[dist+q].I = f4;   f4 = c1[q].R;         f3 = f3 * w[1].R;       \
        f4 = f4 + c1[dist+q].R;                                            \
        f6 = f6 - f5;      f5 = c1[q].I;         f5 = f5 + c1[dist+q].I;   \
        f2 = f2 + f3;      c1[q].R = f4;         c1[dist+q].R = f6;        \
                           c1[q].I = f5;         c1[dist+q].I = f2;
 
#define RISC_DIT(q)                                                        \
      f12 = w[0].R;   f13 = w[0].I;                                        \
      f14 = w[1].R;   f15 = w[1].I;                                        \
      for(i=0;i<q;i++) {                                                   \
         f0 = c0[dist+i].R;  f1 = c0[dist+i].I;                            \
         f2 = c1[dist+i].R;  f3 = c1[dist+i].I;                            \
         f4 = f0 * f12;    f5 = f1 * f13;                                  \
         f8 = f2 * f14;    f9 = f3 * f15;                                  \
         f6 = f0 * f13;    f7 = f1 * f12;  f4 = f4 + f5;                   \
         f10 = f2 * f15;   f11 = f3 * f14; f8 = f8 + f9;                   \
         f6 = f7 - f6;     f0 = c0[i].R;  f10 = f11 - f10; f1 = c0[i].I;   \
         f5 = f4 + f8;     f2 = c1[i].R;  f9 = f4 - f8;    f3 = c1[i].I;   \
         f7 = f6 + f10;    f11 = f6 - f10;                                 \
                                                                           \
         f4 = f0 + f5;   f8 = f0 - f5;                                     \
         f6 = f1 + f7;   f10 = f1 - f7;  c0[i].R = f4;                     \
         f0 = f2 - f11;  c0[dist+i].R = f8;                                \
         f5 = f2 + f11;  c0[i].I = f6;                                     \
         f1 = f3 + f9;   c0[dist+i].I = f10;                               \
         f7 = f3 - f9;   c1[i].R = f0;  c1[dist+i].R = f5;                 \
                         c1[i].I = f1;  c1[dist+i].I = f7;                 \
      }
 
#define PENTIUM_DIT(q)                                                 \
      f0 = c0[dist+q].R;  f1 = w[0].R;  f0 = f0 * f1;                  \
      f1 = c0[dist+q].I;  f2 = w[0].I;  f1 = f1 * f2;                  \
      f2 = c0[dist+q].R;  f3 = w[0].I;  f2 = f2 * f3;                  \
      f3 = c0[dist+q].I;  f4 = w[0].R;  f3 = f3 * f4;                  \
      f4 = c1[dist+q].R;  f5 = w[1].R;  f4 = f4 * f5;                  \
      f5 = c1[dist+q].I;  f6 = w[1].I;  f5 = f5 * f6;                  \
      f6 = c1[dist+q].R;  f7 = w[1].I;  f6 = f6 * f7;                  \
      f7 = c1[dist+q].I;  f7 = f7 * w[1].R;                            \
      f0 = f0 + f1;     f2 = f3 - f2;   f4 = f4 + f5;   f6 = f7 - f6;  \
      f1 = f0;          f1 = f1 - f4;   f0 = f0 + f4;                  \
      f3 = f2;          f3 = f3 - f6;   f2 = f2 + f6;                  \
      f4 = c0[q].R;     f4 = f4 + f0;                                  \
      f5 = c0[q].R;     f0 = f5 - f0;                                  \
      f5 = c0[q].I;     f5 = f5 + f2;  c0[q].R = f4;                   \
      f6 = c0[q].I;     f2 = f6 - f2;  c0[dist+q].R = f0;              \
      f4 = c1[q].R;     f4 = f4 - f3;  c0[q].I = f5;                   \
      f0 = c1[q].R;     f3 = f0 + f3;  c0[dist+q].I = f2;              \
      f5 = c1[q].I;     f5 = f5 + f1;  c1[q].R = f4;                   \
      f2 = c1[q].I;     f1 = f2 - f1;  c1[dist+q].R = f3;              \
                                       c1[q].I = f5;                   \
                                       c1[dist+q].I = f1;
 
/*---------------------------------------------------------------------*/

void single_fft2( cplx *c0, cplx *c1 ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f6,f7;
#else
   double f0,f1,f2,f3,f4;
#endif
 
#ifdef RISC
 
   f0 = c0[0].R;   f1 = c1[0].R;  f2=f0+f1; f6=f0-f1;  
   f0 = c0[0].I;   f1 = c1[0].I;  f3=f0+f1; f7=f0-f1;
   c0[0].R = f2;  c0[0].I = f3;   c1[0].R = f6;  c1[0].I = f7;

#else
 
   f0 = c0[0].R;   f1 = c1[0].R;   f0 = f0 + f1;
   f1 = c0[0].R;   f2 = c1[0].R;   f1 = f1 - f2;
   f2 = c0[0].I;   f3 = c1[0].I;   f2 = f2 + f3;
   f3 = c0[0].I;   f4 = c1[0].I;   f3 = f3 - f4;
     
   c0[0].R = f0;   c1[0].R = f1;
   c0[0].I = f2;   c1[0].I = f3;

#endif
 
}
 
/*---------------------------------------------------------------------*/

void vector_fft2( cplx *c0, cplx *c1 ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9;
#else
   double f0,f1,f2,f3,f4;
#endif
 
#ifdef RISC
 
   f0 = c0[0].R;   f1 = c1[0].R;  f2=f0+f1; f6=f0-f1;  
   f0 = c0[0].I;   f1 = c1[0].I;  f3=f0+f1; f7=f0-f1;
   f0 = c0[1].R;   f1 = c1[1].R;  f4=f0+f1; f8=f0-f1;
   f0 = c0[1].I;   f1 = c1[1].I;  f5=f0+f1; f9=f0-f1;
   c0[0].R = f2;  c0[0].I = f3;   c1[0].R = f6;  c1[0].I = f7;
   c0[1].R = f4;  c0[1].I = f5;   c1[1].R = f8;  c1[1].I = f9;
 
   f0 = c0[2].R;   f1 = c1[2].R;  f2=f0+f1; f6=f0-f1;  
   f0 = c0[2].I;   f1 = c1[2].I;  f3=f0+f1; f7=f0-f1;
   f0 = c0[3].R;   f1 = c1[3].R;  f4=f0+f1; f8=f0-f1;
   f0 = c0[3].I;   f1 = c1[3].I;  f5=f0+f1; f9=f0-f1;
   c0[2].R = f2;  c0[2].I = f3;   c1[2].R = f6;  c1[2].I = f7;
   c0[3].R = f4;  c0[3].I = f5;   c1[3].R = f8;  c1[3].I = f9;
 
#else
 
   f0 = c0[0].R;   f1 = c1[0].R;   f0 = f0 + f1;
   f1 = c0[0].R;   f2 = c1[0].R;   f1 = f1 - f2;
   f2 = c0[0].I;   f3 = c1[0].I;   f2 = f2 + f3;
   f3 = c0[0].I;   f4 = c1[0].I;   f3 = f3 - f4;
   c0[0].R = f0;   f0 = c0[1].R;   f0 += c1[1].R;
   c1[0].R = f1;   f1 = c0[1].R;   f1 -= c1[1].R;
   c0[0].I = f2;   f2 = c0[1].I;   f2 += c1[1].I;
   c1[0].I = f3;   f3 = c0[1].I;   f3 -= c1[1].I;
   c0[1].R = f0;   f0 = c0[2].R;   f0 += c1[2].R;
   c1[1].R = f1;   f1 = c0[2].R;   f1 -= c1[2].R;
   c0[1].I = f2;   f2 = c0[2].I;   f2 += c1[2].I;
   c1[1].I = f3;   f3 = c0[2].I;   f3 -= c1[2].I;
   c0[2].R = f0;   f0 = c0[3].R;   f0 += c1[3].R;
   c1[2].R = f1;   f1 = c0[3].R;   f1 -= c1[3].R;
   c0[2].I = f2;   f2 = c0[3].I;   f2 += c1[3].I;   c0[3].R = f0;
   c1[2].I = f3;   f3 = c0[3].I;   f3 -= c1[3].I;   c1[3].R = f1;
                                                    c0[3].I = f2;
                                                    c1[3].I = f3;
#endif
 
}
 
/*---------------------------------------------------------------------*/

void single_splitrad_DIF( cplx *c, cplx *w, long size ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6;
#endif

   long dist = size/4;
   cplx *c0 = c,
        *c1 = c + dist;
 
   dist = 2*dist;
   size = size/4;
 
#ifdef RISC
 
   do {
      RISC_DIF(1);
        
      c0++;  c1++;
      w+=2;   size--;      
   } while(size);
 
#else
 
   do {
        PENTIUM_DIF(0)
 
        c0++;  c1++;
        w+=2;   size--;     
   } while(size);
 
#endif
 
}
 
/*---------------------------------------------------------------------*/

void vector_splitrad_DIF( cplx *c, cplx *w, long size ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6;
#endif

   long dist = size/4 * 4;
   cplx *c0 = c,
        *c1 = c + dist;
 
   dist = 2*dist;
   size = size/4;
 
#ifdef RISC
 
   do {
      RISC_DIF(4);
        
      c0+=4;  c1+=4;
      w+=2;   size--;      
   } while(size);
 
#else
 
   do {
        PENTIUM_DIF(0)
        PENTIUM_DIF(1)
        PENTIUM_DIF(2)
        PENTIUM_DIF(3)
 
        c0+=4;  c1+=4;
        w+=2;   size--;     
   } while(size);
 
#endif
 
}
 
/*---------------------------------------------------------------------*/

void single_splitrad_DIT( cplx *c, cplx *w, long size ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif

   long dist = size/4;
   cplx *c0 = c,
        *c1 = c + dist;
 
   dist = 2*dist;
   size = size/4;
 
#ifdef RISC
 
   do {
      RISC_DIT(1)
 
      c0++;   c1++;
      w+=2;   size--;      
   } while(size);
 
#else
 
   do {
      PENTIUM_DIT(0)
 
      c0++;   c1++;
      w+=2;   size--;      
   } while(size);
 
#endif
 
}
 
/*---------------------------------------------------------------------*/

void vector_splitrad_DIT( cplx *c, cplx *w, long size ) {
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
   long   i;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif

   long dist = size/4 * 4;
   cplx *c0 = c,
        *c1 = c + dist;
 
   dist = 2*dist;
   size = size/4;
 
#ifdef RISC
 
   do {
      RISC_DIT(4)
 
      c0+=4;  c1+=4;
      w+=2;   size--;      
   } while(size);
 
#else
 
   do {
      PENTIUM_DIT(0)
      PENTIUM_DIT(1)
      PENTIUM_DIT(2)
      PENTIUM_DIT(3)
 
      c0+=4;  c1+=4;
      w+=2;   size--;      
   } while(size);
 
#endif
 
}
 
/*---------------------------------------------------------------------*/

void base1(cplx *c) {
 
        /* base of the recursion for the first pass:
         *
         *  x[0]     x1       (x1 + z1)^2
         *  x[1] --> z1 -->   (x1 - z1)^2
         *  x[2]     x2     (x1*x2 - z1*z2)
         *  x[3]     z2     (x1*z2 - z1*x2)
         */
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif
 
#ifdef RISC
 
   f0 = c[0].R;      f1 = c[0].I;  
   f2 = c[1].R;      f3 = c[1].I;    f8 = f0 + f2;  f9 = f0 - f2;
   f4 = c[2].R;      f5 = c[2].I;    f10 = f1 + f3; f11 = f1 - f3;
   f6 = c[3].R;      f7 = c[3].I;    f12 = f8 * f10; f13 = f8 * f8;
   f14 = f10 * f10; f12 = f12 + f12; f8 = f9 * f11; f13 = f13 - f14;
   f10 = f9 * f9;   c[0].I = f12;    c[0].R = f13;  f11 = f11 * f11;
   f8 = f8 + f8;    f12 = f0 * f4;   f13 = f1 * f5;
   f14 = f0 * f5;   f10 = f10 - f11;  c[1].I = f8;  f15 = f1 * f4;
   f8 = f2 * f6;    f9 = f3 * f7;     f12 = f12 - f13;  c[1].R = f10;
   f10 = f2 * f7;   f14 = f14 + f15;  f11 = f3 * f6;    f8 = f8 - f9;
   f13 = f1 * f7;   f15 = f1 * f6;    f10 = f10 + f11;  f8 = f12 - f8;
   f9 = f14 - f10;   c[2].R = f8;      c[2].I = f9;     f14 = f0 * f7;   
   f12 = f0 * f6;   f8 = f2 * f4;     f9 = f3 * f5;     f14 = f14 + f15;
   f10 = f2 * f5;   f12 = f12 - f13;  f11 = f3 * f4;    f8 = f8 - f9;
   f10 = f10 + f11; f8 = f12 - f8;    f10 = f14 - f10;
   c[3].R = f8;     c[3].I = f10;
 
#else
 
   f0 = c[0].R;  f1 = c[2].R;  f0 = f0 * f1;
   f1 = c[0].I;  f2 = c[2].I;  f1 = f1 * f2;
   f2 = c[0].R;  f3 = c[2].I;  f2 = f2 * f3;
   f3 = c[0].I;  f4 = c[2].R;  f3 = f3 * f4;
   f4 = c[1].R;  f5 = c[3].R;  f4 = f4 * f5;
   f5 = c[1].I;  f6 = c[3].I;  f5 = f5 * f6;
   f6 = c[1].R;  f7 = c[3].I;  f6 = f6 * f7;
   f7 = c[1].I;  f7 = f7 * c[3].R;
   f0 = f0 - f1; f4 = f4 - f5; f2 = f2 + f3; f6 = f6 + f7;
   f0 = f0 - f4; f2 = f2 - f6;
   f4 = c[2].R;  f5 = c[1].R;  f4 = f4 * f5;
   f5 = c[2].I;  f6 = c[1].I;  f5 = f5 * f6;
   f6 = c[2].R;  f7 = c[1].I;  f6 = f6 * f7;
   f7 = c[2].I;  f7 = f7 * c[1].R;  c[2].R = f0; c[2].I = f2; 
   f0 = c[0].R;  f1 = c[3].R;  f0 = f0 * f1;
   f1 = c[0].I;  f2 = c[3].I;  f1 = f1 * f2;
   f2 = c[0].R;  f3 = c[3].I;  f2 = f2 * f3;
   f3 = c[0].I;  f3 = f3 * c[3].R;
   f4 = f4 - f5; f6 = f6 + f7; f0 = f0 - f1; f2 = f2 + f3;
   f0 = f0 - f4; f1 = c[0].R;  f3 = c[1].R;  f1 = f1 + f3;
   f2 = f2 - f6; c[3].R = f0;  f0 = c[0].I;  f3 = c[1].I;  f0 = f0 + f3;
   c[3].I = f2;  f2 = c[0].I;  f3 = c[1].I;  f2 = f2 - f3;
   f3 = c[0].R;  f4 = c[1].R;  f3 = f3 - f4;
   f4 = f0;      f4 = f4 * f1; f1 = f1 * f1; f0 = f0 * f0;
   f5 = f2;      f5 = f5 * f3; f4 = f4 + f4; f2 = f2 * f2;
   f1 = f1 - f0; f3 = f3 * f3; f5 = f5 + f5; c[0].I = f4;
   f3 = f3 - f2; c[0].R = f1;  c[1].I = f5;  c[1].R = f3;
 
#endif
 
}
 
/*---------------------------------------------------------------------*/

void twist( cplx *c, long size, long vec ) {
 
   /* applies "size" forward dwt factors, each to "vec"
      complex's at a time.             */
 
   cplx *d = dwt;
   long i;
   double f0,f1,f2,f3,f4,f5;
 
   do {
      f4 = d[0].R;   f5 = d[0].I;
      for(i=0; i<vec; i++,c++) {
         f0 = c[0].R;  f0 *= f4;
         f1 = c[0].I;  f1 *= f5;
         f2 = c[0].R;  f2 *= f5;
         f3 = c[0].I;  f3 *= f4;
         f0 = f0 - f1; f2 = f2 + f3;
         c[0].R = f0;  c[0].I = f2;
      }
      d++; size--;
   } while(size);
}
 
/*---------------------------------------------------------------------*/

void untwist( cplx *c, long size, long vec ) {
 
   /* same as "twist", but conjugates the dwt factors
      before multiplying */
 
   cplx *d = dwt;
   long i;
   double f0,f1,f2,f3,f4,f5;
 
   do {
      f4 = d[0].R;   f5 = d[0].I;
      for(i=0; i<vec; i++,c++) {
         f0 = c[0].R;  f0 *= f4;
         f1 = c[0].I;  f1 *= f5;
         f2 = c[0].R;  f2 *= f5;
         f3 = c[0].I;  f3 *= f4;
         f0 = f0 + f1; f2 = f3 - f2;
         c[0].R = f0;  c[0].I = f2;
      }
      d++; size--;
   } while(size);
}
 
/*---------------------------------------------------------------------*/

void release_carries( cplx *c, long size, double scale, long vec ) {
 
   /* release size carries for vec arrays of cplx data,
      multiplying by scale in the process. The real and
      imaginary carries are released separately. To
      fully finish releasing carries use last_carry() */
 
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,
          f9,f10;
   long i;
 
   f2 = magic[0];            /* magic rounding constant */
   f3 = magic[1];            /* same (to fool optimizing compilers) */
   f4 = magic[2];            /* 1/BASE */
   f5 = magic[3];            /* BASE */
   f10 = scale;
 
   do {
      for(i=0; i<vec; i++,c++) {
         f0 = carries[i].R; f1 = carries[i].I;
         f6 = c[0].R;   f7 = c[0].I;
         f6 = f6 * f10; f7 = f7 * f10;  /* scale the FFT result */
         f6 += f0;      f7 += f1;       /* add previous carries */
         f6 += f2;      f7 += f2;       /* round the result */
         f6 -= f3;      f7 -= f3;
         f8 = f6 * f4;  f9 = f7 * f4;   /* divide by BASE */
         f8 += f2;      f9 += f2;       /* round the result */
         f8 -= f3;      f9 -= f3;
         carries[i].R = f8;
         carries[i].I = f9;             /* new carries */
         f8 *= f5;      f9 *= f5;
         f6 -= f8;      f7 -= f9;       /* compute the remainder */
         c[0].R = f6;   c[0].I = f7;
      }
      size--;
   }while(size);
 
}
 
/*---------------------------------------------------------------------*/

void last_carry( cplx *c, long vec ) {
 
   /* after calling release_carries(), call this on the
      beginning of the array to perform a Fermat mod on the
      carry and propagate the carry again until 0. If
      10 elements propagated aren't enough to zero out all
      the carries this routine aborts the program */
 
   double f0, f1;
   long i, length = (Q/32<10)?(Q/32):(10);
 
   for(i=0; i<vec; i++) {
      f0 = carries[i].R;
      f1 = carries[i].I;
      carries[i].R = -f1;
      carries[i].I = f0;
   }
 
   release_carries(c, length, 1.0, vec );
 
   for(i=0; i<vec; i++) {
      if( carries[i].R || carries[i].I ) {
         printf("\nerror: carry not zero!\n");
         exit(-1);
      }
   }
}
