/**************************************************************

    ellmul0.c
   
    Separate source for Windows compilation, to get around
    Visual C++ 5.0 compiler bug.

    Change history:

    28 Jul 99   ARP - Added conditional variable declarations
    28 Jul 99   REC - Creation.

    c. 1999 Perfectly Scientific, Inc.
    All Rights Reserved.

**************************************************************/

typedef struct {
   double R, I;
} cplx;

/*----------------- global variables ------------------------*/

extern cplx *cacheptr, *temptr;

/*--------------- local prototypes --------------------------*/

void base2(cplx *c);
void base3(cplx *c);

/*-----------------------------------------------------------*/

void base2(cplx *c) {
 
        /* base of the recursion for the second pass:
         *     An = cacheptr[0];
         *     Ad = cacheptr[1];
         *
         *  x0            x0 * x1
         *  x1 -->   Ad * x1 + An * (x0-x1)
         *  x2            x2 * x2
         *  x3            x3 * x3 
         */
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif
 
#ifdef RISC
 
   f0 = c[3].R;  f1 = c[3].I;
   f2 = c[2].R;  f3 = c[2].I;  f8 = f0 * f1;  f9 = f0 * f0;
   f4 = c[1].R;  f5 = c[1].I;  f10 = f1 * f1; f8 = f8 + f8;
   f6 = c[0].R;  f7 = c[0].I;  f11 = f2 * f3;  f12 = f2 * f2;
   f9 = f9 - f10;  f13 = f3 * f3;  f11 = f11 + f11;
   c[3].I = f8;    f14 = f6 - f4;  c[3].R = f9;  f15 = f7 - f5;
   c[2].I = f11;   f12 = f12 - f13; c[2].R = f12;
   f0 = cacheptr[1].R;  f1 = cacheptr[1].I;  f8 = f4 * f6;
   f2 = cacheptr[0].R;  f3 = cacheptr[0].I;  f9 = f5 * f7;
   f10 = f4 * f7;  f11 = f5 * f6; f12 = f2 * f14;  f13 = f3 * f15;
   f8 = f8 - f9;   f6 = f0 * f4; f10 = f10 + f11; f7 = f1 * f5;
   c[0].R = f8;    f14 = f3 * f14;  c[0].I = f10;  f15 = f2 * f15;
   f12 = f12 - f13;  f0 = f0 * f5;  f6 = f6 - f7;  f1 = f1 * f4;
   f14 = f14 + f15;  f0 = f0 + f1;  f12 = f12 + f6;  c[1].R = f12;
   f14 = f14 + f0;   c[1].I = f14;
 
   cacheptr += 4;
 
#else
 
   f0 = c[2].R;  f0 = f0 * f0;
   f1 = c[2].I;  f1 = f1 * f1;
   f2 = c[3].R;  f2 = f2 * f2;
   f3 = c[3].I;  f3 = f3 * f3;
   f4 = c[2].R;  f5 = c[2].I;   f4 = f4 * f5;
   f5 = c[3].R;  f6 = c[3].I;   f5 = f5 * f6;
   f0 = f0 - f1;  f2 = f2 - f3;  f4 = f4 + f4;  f5 = f5 + f5;
   c[2].R = f0;   c[3].R = f2;   c[2].I = f4;   c[3].I = f5;
   f0 = c[0].R;  f1 = c[1].R;   f0 = f0 - f1;
   f1 = c[0].I;  f2 = c[1].I;   f1 = f1 - f2;
   f2 = cacheptr[0].R;  f2 = f2 * f0;
   f3 = cacheptr[0].I;  f3 = f3 * f1;
   f4 = cacheptr[0].I;  f0 = f4 * f0;
   f4 = cacheptr[0].R;  f1 = f4 * f1;
   f4 = cacheptr[1].R;  f5 = c[1].R;  f4 = f4 * f5;
   f5 = cacheptr[1].I;  f6 = c[1].I;  f5 = f5 * f6;
   f6 = cacheptr[1].R;  f7 = c[1].I;  f6 = f6 * f7;
   f7 = cacheptr[1].I;  f7 = f7 * c[1].R;
   f2 = f2 - f3;   f0 = f0 + f1;   f4 = f4 - f5;   f6 = f6 + f7;
   f2 = f2 + f4;   f1 = c[0].R;    f3 = c[1].R;    f1 = f1 * f3;
   f0 = f0 + f6;   f3 = c[0].I;    f5 = c[1].I;    f3 = f3 * f5;
   f5 = c[0].R;    f7 = c[1].I;    f5 = f5 * f7;
   f7 = c[0].I;    f7 = f7 * c[1].R;  c[1].R = f2;  f1 = f1 - f3;
   c[1].I = f0;    f5 = f5 + f7;   c[0].R = f1;     c[0].I = f5;
 
   cacheptr += 4;
 
#endif
 
}
 
/*-----------------------------------------------------------*/

void base3(cplx *c) {
 
        /* base of the recursion for the third pass:
         *
         *  x[0]     x[0] * cacheptr[1]
         *  x[1] --> x[1] * temptr[0]
         *  x[2]     x[2] * cacheptr[2]
         *  x[3]     x[3] * cacheptr[3]
         */
 
#ifdef RISC
   double f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14,f15;
#else
   double f0,f1,f2,f3,f4,f5,f6,f7;
#endif
 
#ifdef RISC
 
   f0 = c[0].R;         f1 = c[0].I;
   f2 = cacheptr[1].R;  f3 = cacheptr[1].I;
   f4 = c[1].R;         f5 = c[1].I;        f8 = f0 * f2;  f9 = f1 * f3;
   f6 = temptr[0].R;    f7 = temptr[0].I;   f10 = f0 * f3;  f11 = f1 * f2;
   f8 = f8 - f9;     f12 = f4 * f6;   f13 = f5 * f7;   f14 = f4 * f7;
   f10 = f10 + f11;  f12 = f12 - f13;  f15 = f5 * f6;  c[0].R = f8;
   c[0].I = f10;     c[1].R = f12;     f14 = f14 + f15;  c[1].I = f14;
   f0 = c[2].R;         f1 = c[2].I;
   f2 = cacheptr[2].R;  f3 = cacheptr[2].I;
   f4 = c[3].R;         f5 = c[3].I;        f8 = f0 * f2;  f9 = f1 * f3;
   f6 = cacheptr[3].R;  f7 = cacheptr[3].I; f10 = f0 * f3;  f11 = f1 * f2;
   f8 = f8 - f9;     f12 = f4 * f6;   f13 = f5 * f7;   f14 = f4 * f7;
   f10 = f10 + f11;  f12 = f12 - f13;  f15 = f5 * f6;  c[2].R = f8;
   c[2].I = f10;     c[3].R = f12;     f14 = f14 + f15;  c[3].I = f14;
 
   cacheptr += 4; temptr++;
 
#else
 
   f0 = c[0].R;  f1 = cacheptr[1].R;  f0 = f0 * f1;
   f1 = c[0].I;  f2 = cacheptr[1].I;  f1 = f1 * f2;
   f2 = c[0].R;  f3 = cacheptr[1].I;  f2 = f2 * f3;
   f3 = c[0].I;  f4 = cacheptr[1].R;  f3 = f3 * f4;
   f4 = c[1].R;  f5 = temptr[0].R;  f4 = f4 * f5;
   f5 = c[1].I;  f6 = temptr[0].I;  f5 = f5 * f6;
   f6 = c[1].R;  f7 = temptr[0].I;  f6 = f6 * f7;
   f7 = c[1].I;  f7 = f7 * temptr[0].R;
   f0 = f0 - f1; f4 = f4 - f5; f2 = f2 + f3; f6 = f6 + f7;
   c[0].R = f0;  c[1].R = f4;   c[0].I = f2;  c[1].I = f6;
   f0 = c[2].R;  f1 = cacheptr[2].R;  f0 = f0 * f1;
   f1 = c[2].I;  f2 = cacheptr[2].I;  f1 = f1 * f2;
   f2 = c[2].R;  f3 = cacheptr[2].I;  f2 = f2 * f3;
   f3 = c[2].I;  f4 = cacheptr[2].R;  f3 = f3 * f4;
   f4 = c[3].R;  f5 = cacheptr[3].R;  f4 = f4 * f5;
   f5 = c[3].I;  f6 = cacheptr[3].I;  f5 = f5 * f6;
   f6 = c[3].R;  f7 = cacheptr[3].I;  f6 = f6 * f7;
   f7 = c[3].I;  f7 = f7 * cacheptr[3].R;
   f0 = f0 - f1; f4 = f4 - f5; f2 = f2 + f3; f6 = f6 + f7;
   c[2].R = f0;  c[3].R = f4;   c[2].I = f2;  c[3].I = f6;
 
   cacheptr += 4; temptr++;
 
#endif
 
}
