/**************************************************************
    fermat.c

    Automated attack program for large Fermat numbers.

    Updates:
        05 Dec 99   REC - Installed J. Papadopoulos' fix pr[i]*pr[i]
        27 Jul 99   REC - Incorporated some fine tunings by Powell
        21 Jul 99   REC - More safeguards on allowed B limits.
        21 Jul 99   REC - Network version with safeguards & tests.
        18 Jul 99   REC - Fixed special_mulg() bug
        17 Apr 99   REC - Comment: McIntosh/Tardif factor of F_18
        05 Mar 99   REC - Memory improvements
        25 Feb 99   REC - Reduced variable storage
        23 Feb 99   REC - Various improvements a la McIntosh, Powell
        24 Jul 97   RDW - added new note text
        20 May 97   RDW - fixed Win32 compiler warnings
        18 May 97   REC - changed s_mod()
        26 Apr 97	RDW - fixed tabs and unix EOL
        24 Apr 97	REC
        19 Apr 97   RDW

    c. 1997 Perfectly Scientific, Inc.
    All Rights Reserved.

    To compile (note some platforms admit of various cc or gcc
    options):
    
    Windows: 
 
    % cl ellmul0.c -c
    % cl /Ox /G5 /Op giants.c fermat.c ellmul.c ellmul0.obj -o fermat.exe
    _________________________________
    Linux:

    % gcc -O3 fermat.c giants.c ellmul.c ellmul0.c -lm
-o fermat
    % gcc -O3 fermat.c giants.c ellmul.c ellmul0.c -m486 -fomit-frame-pointer -fno-defer-pop -lm

   _________________________________
    OpenStep:

    % cc -O fermat.c giants.c ellmul.c ellmul0.c -lm -o fermat
    _________________________________
    SGI:

    % cc -O fermat.c giants.c ellmul.c ellmul0.c -lm -DRISC -o fermat 
    _________________________________
    Sun:

    % gcc -O3 fermat.c giants.c ellmul.c ellmul0.c -lm -DRISC -o fermat
    % gcc -O3 -funroll-loops -msupersparc -ffast-math fermat.c giants.c ellmul.c ellmul0.c -lm
-DRISC -o fermat
    % cc -fast -native -xO5 -dalign -xspace fermat.c giants.c ellmul.c ellmul0.c -lm -DRISC -o fermat
    _________________________________
    DEC Alpha:

    % cc -O fermat.c giants.c ellmul.c ellmul0.c -lm -o fermat -DRISC
    _________________________________

    Discussion:

    The basic algorithm is the elliptic curve method (ECM) of 
    H. Lenstra, with
    R. Brent's seed parameterization.  The discrete weighted
    transform (DWT) of R. Crandall and B. Fagin is used
    for the big multiplies.  Various ECM enhancements due to
    R. Brent, R. Crandall, P. Montgomery, G. Woltman, 
    and P. Zimmerman are also used.

    Credit is due also to: K. Dilchar, R. McIntosh for valuable testing
    and optimizations; to A. Powell for substantial improvements
    to the companion large-integer library 'giants.c', and to
    J. P. Buhler for various theoretical observations over the years.
    The present 'network' version of fermat.c has been enhanced
    dramatically by J. Papadopoulos, especially in the form of the
    companion library 'ellmul.c' which implements many of the
    aforementioned modern enhancements.

    Usage:

    To attack F_n = 2^(2^n) + 1 with recommended internal paremeters:

    % fermat n
 
    To test performance/integrity on your machine:

    % fermat n -t

    (which runs known elliptic curves, with certain expected results).
    To force a stage limit of B = 1000000, say:

    % fermat n -b 1000000 

    where the algorithm is appropriate for n one of:
    7 through 22, inclusive. 
    A random ECM seed is determined, not only
    for every new runtime, but subsequent seeds during any one
    long run are also random.

    The program outputs a file whose progress, in choosing
    random ECM seeds, and carrying out stage one and stage
    two of ECM, is signified by dots ".". These dots emanate
    slowly.  Some machines will take several days to complete
    an ECM seed.  One can look at the output file to check on
    the printing of dots.
    If a factor is found it (together with parameter values)
    is prefixed with the character '!'
    (for stage-one discovery), or '!!' for stage-two discovery.
    Thus, by searching for '!', for example on Unix systems via

    % grep '\!' file(s)

    one can hope for success in the form of lines prefixed with one
    or more '!' characters.
    The program is working fine if it produces the known
    F_n factors (as of Jul 97):

    Notes on known factors of F_n:

    -- 1991-2, the second and third listed factors of F_13, below,
       were found by R. Crandall using this very algorithm.
       Later, Brent found the largest known factor of F_13 using
       his variant (after which ours ia modelled).

    -- 4 Jan 97, Crandall & Dilcher found the new
       (second) factor of F_16 using this very algorithm.

    -- 3 Jul 97, Crandall & Van Halewyn found the new
       (33-digit) factor of F_15 using this very algorithm.

    -- Every F_n thgrough n = 11 has been completely
       factored.  For example, R. Brent demolished F_10 in 1996,
       and F_9 was resolved via the NFS also in the 1990's.

    -- McIntosh and Tardif found the 23-digit factor of F_18
       on 16 Apr 1999).

    -- A (proven) composite number is called GENUINE if no nontrivial factor
       be known.  All the listed 'composites' below are genuine.
       F_20 was proved composite by Young & Buell.  F_22 and various of the
       other 'composite' factors below were proven so by Crandall, Doenias, Norrie,
       Young.

    -- A more complete summary of Fermat-number status is at W. Keller's page:
       http://vamri.xray.ufl.edu/proths/fermat.html

    F0-F4: prime
    F5-F11: completely factored
    F12 = 114689 * 26017793 * 63766529 * 190274191361 *  	    	
          1256132134125569 * composite
    F13 = 2710954639361 * 2663848877152141313 * 3603109844542291969 *
          319546020820551643220672513 * composite
    F14 = composite
    F15 = 1214251009 * 2327042503868417 *
              168768817029516972383024127016961 * composite
    F_16 = 825753601 * 188981757975021318420037633 * composite
    F_17 = 31065037602817 * composite
    F_18 = 13631489 * 81274690703860512587777 * composite
    F_19 = 70525124609 * 646730219521 * composite
    F_20 = composite
    F_21 = 4485296422913 * composite
    F_22 = composite

    When running this program, beware of the simultaneous discovery of combinations,
    such as *both* of the known F_19
    factors multiplied together...that might appear on first glance
    as a victory!
 **************************************************************/


/* include files */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#ifdef _WIN32
#include <process.h>
#endif
#include "giants.h"


/* definitions */

#define s_modg(N, g) fer_mod(Q, g, N)
#define D_MAX 600
#define NUM_PRIMES 6542 	/* PrimePi[2^16]. */

/* compiler options */

#ifdef _WIN32
#pragma warning( disable : 4127 4706 ) /* disable conditional is constant warning */
#endif


/* global variables */

double 		*xf, *yf;
int 		D, pr[NUM_PRIMES], pjump[NUM_PRIMES];
unsigned char parr[D_MAX];
giant		xr = NULL, xs = NULL, zs = NULL, zr = NULL,
			xorg = NULL, zorg = NULL, t1 = NULL, t2 = NULL,
			t3 = NULL, N = NULL, gg = NULL, An = NULL, Ad = NULL;
giant		xb[D_MAX+2], zb[D_MAX+2], xzb[D_MAX+2];
int			Q;

extern giant    cur_den, cur_recip, scratch;


/* function prototypes */

int 		least_2(int n);
int             psi_rand(void);

/* function prototypes from ellmul.c */

void            init_ellmul(long size); 
void            ell_mul( giant xx, giant zz, int n, giant An,
                          giant Ad, giant N );
void            special_mulg( giant y, giant x );
void            special_squareg( giant y );



/**************************************************************
 *
 *	Functions
 *
 **************************************************************/

void
set_random_seed(
	void
)
{
	/* Start the random number generator at a new position. */
	time_t	tp;

	time(&tp);
  	srand((int)tp + (int)getpid());
}


int			
psi_rand(
	void
)
{
	unsigned short	hi;
	unsigned short	low;
	time_t			tp;
	int				result;

	time(&tp);
	low = (unsigned short)rand();
	hi = (unsigned short)rand();
	result = ((hi << 16) | low) ^ ((int)tp);

	return (result & 0x7fffffff);
}


int
isprime(
	unsigned int	odd
)
{
	int 			j;
	unsigned int 	p;

	for (j=1; ; j++)
	{
		 p = pr[j];
		 if (p*p > odd)
			return(1);
		 if (odd % p == 0)
			return(0);
	}
}



int
primeq(
	int				p
)
{
	register int	j=3;

	if ((p&1)==0)
		return((p==2)?1:0);
	if (j>=p)
		return(1);
	while ((p%j)!=0)
	{
		j+=2;
		if (j*j>p)
			return(1);
	}
	return(0);
}


int
power_2(
	int				s
)
{
	unsigned int 	c = (unsigned int)(1<<31);
	int 			p = 31;

	while (c)
	{
		if(c==(unsigned int)s)
			return(p);
		--p;
		c >>= 1;
	}
	return(-1);
}


void
ell_even(
	giant 	x1,
	giant  	z1,
	giant  	x2,
	giant  	z2,
	giant  	An,
	giant  	Ad,
	giant  	N
)
{
	gtog(x1, t1);
	addg(z1, t1);
	special_squareg(t1);
	s_modg(N, t1);
	gtog(x1, t2);
	subg(z1, t2);
	special_squareg(t2);
	s_modg(N, t2);
	gtog(t1, t3);
	subg(t2, t3);
	gtog(t2, x2);
	special_mulg(t1, x2);
	gshiftleft(2, x2);
	s_modg(N, x2);
	special_mulg(Ad, x2);
	s_modg(N, x2);
	special_mulg(Ad, t2);
	gshiftleft(2, t2);
	s_modg(N, t2);
	gtog(t3, t1);
	special_mulg(An, t1);
	s_modg(N, t1);
	addg(t1, t2);
	special_mulg(t3, t2);
	s_modg(N, t2);
	gtog(t2,z2);
}



void
ell_odd(
	giant 	x1,
	giant  	z1,
	giant 	x2,
	giant  	z2,
	giant 	xor,
	giant 	zor,
	giant  	N
)
{
	gtog(x1, t1);
	subg(z1, t1);
	gtog(x2, t2);
	addg(z2, t2);
	special_mulg(t1, t2);
	s_modg(N, t2);
	gtog(x1, t1);
	addg(z1, t1);
	gtog(x2, t3);
	subg(z2, t3);
	special_mulg(t3, t1);
	s_modg(N, t1);
	gtog(t2, x2);
	addg(t1, x2);
	special_squareg(x2);
	s_modg(N, x2);
	gtog(t2, z2);
	subg(t1, z2);
	special_squareg(z2);
	s_modg(N, z2);
	special_mulg(zor, x2);
	s_modg(N, x2);
	special_mulg(xor, z2);
	s_modg(N, z2);
}


/* From R. P. Brent, priv. comm. 1996:
Let s > 5 be a pseudo-random seed (called $\sigma$ in the Tech. Report),

	u/v = (s^2 - 5)/(4s)

Then starting point is (x_1, y_1) where

	x_1 = (u/v)^3
and
	a = (v-u)^3(3u+v)/(4u^3 v) - 2
*/



void
choose12
(
	giant 	x,
	giant  	z,
	int    	k,
	giant 	An,
	giant 	Ad,
	giant 	N
)
{
	itog(k, zs);
	gtog(zs, xs);
	special_squareg(xs);
	itog(5, t2);
	subg(t2, xs);
	s_modg(N, xs);
	addg(zs, zs);
	addg(zs, zs);
	s_modg(N, zs);
	gtog(xs, x);
	special_squareg(x);
	s_modg(N, x);
	special_mulg(xs, x);
	s_modg(N, x);
	gtog(zs, z);
	special_squareg(z);
	s_modg(N, z);
	special_mulg(zs, z);
	s_modg(N, z);

	/* Now for A. */
	gtog(zs, t2);
	subg(xs, t2);
	gtog(t2, t3);
	special_squareg(t2);
	s_modg(N, t2);
	special_mulg(t3, t2);
	s_modg(N, t2);  /* (v-u)^3. */
	gtog(xs, t3);
	addg(t3, t3);
	addg(xs, t3);
	addg(zs, t3);
	s_modg(N, t3);
	special_mulg(t3, t2);
	s_modg(N, t2);  /* (v-u)^3 (3u+v). */
	gtog(zs, t3);
	special_mulg(xs, t3);
	s_modg(N, t3);
	special_squareg(xs);
	s_modg(N, xs);
	special_mulg(xs, t3);
	s_modg(N, t3);
	addg(t3, t3);
	addg(t3, t3);
	s_modg(N, t3);
	gtog(t3, Ad);
	gtog(t2, An);  /* An/Ad is now A + 2. */
}



void
ensure(
	int 	q
)
{
	int 	nsh, j;

	if (!q)
	{
		gread(N,stdin);
		q = bitlen(N) + 1;
	}

        nsh = Q/16 + 32;

	N = newgiant(nsh);

        if (!xr) xr = newgiant(nsh);
        if (!zr) zr = newgiant(nsh);
        if (!xs) xs = newgiant(nsh);
	if (!zs) zs = newgiant(nsh);
	if (!xorg) xorg = newgiant(nsh);
	if (!zorg) zorg = newgiant(nsh);
	if (!t1) t1 = newgiant(nsh);
	if (!t2) t2 = newgiant(nsh);
	if (!t3) t3 = newgiant(nsh);
	if (!gg) gg = newgiant(nsh);
	if (!An) An = newgiant(nsh);
        if (!Ad)  Ad = newgiant(nsh);

	for (j=0;j<D+2;j++)
	{
		xb[j] = newgiant(nsh);
		zb[j] = newgiant(nsh);
		xzb[j] = newgiant(nsh);
	}

        init_ellmul(Q/16);

}



int
least_2(
	int				n
)
{
	/* This function returns least power of two greater than n. */
	register int	i = 1;

	while(i<n)
	{
		i<<=1;
	}
	return(i);
}



void
dot(
	void
)
{
	printf(".");
	fflush(stdout);
}



/**************************************************************
 *
 *	Main Function
 *
 **************************************************************/

main(
	int				argc,
	char 			*argv[]
)
{
	unsigned int	i, iend, j, k, m, pri;
	int 			m1, m2, test_mode;
	int             nshorts, cnt, count, limitbits = 0, pass, npr, stage_flag, seed;
	unsigned long	B, C;

	if (argc < 2)
	{
		fprintf(stdout,"Usage:\nTo attack F_n with recommended internal parameters:\n%% fermat n\n");
		fprintf(stdout,"To test the performance/integrity on your machine:\n%% fermat n -t\n");
		fprintf(stdout,"To force stage-one limit parameter:\n%% fermat n -b <limit>\n");
		exit(0);
	}
    if((argc == 3) && (!strcmp(argv[2], "-t"))) test_mode = 1;
		else test_mode = 0;
    Q = atoi(argv[1]);
	if((Q < 7) || (Q > 22)) {
		fprintf(stderr, "Only those F_n with 7 <= n <= 22 allowed.\n");
		exit(0);
	}
	Q = 1 << Q;
	printf("Attacking 2^%d + 1\n", Q);
	fflush(stdout);
	switch (atoi(argv[1]))
	{
        case 12: B = 20000000L; D = 512; break;
        case 13: B =  8000000L; D = 512; break;
		case 14: B =  2000000L; D = 512; break;
		case 15: B =  1000000L; D = 256; break;
		case 16: B =   400000L; D = 256; break;
		case 17: B =   200000L; D = 256; break;
		case 18: B =   100000L; D = 128; break;
		case 19: B =    40000L; D = 64; break;
		case 20: B =    20000L; D = 32; break;
		case 21: B =    10000L; D = 16; break;
		case 22: B =    10000L; D =  8; break;
		default: B = 10000000L; D = 512; break;
	}

    if((argc == 4) && (!strcmp(argv[2], "-b"))) {
		B = atoi(argv[3]);
        if((B < 100) || (B > 40000000)) {
			fprintf(stderr, "B limit must be between 100, 40000000 inclusive\n");
			exit(0);
		}
	}
	C = 50*B;
	ensure(Q);
	itog(1, N);
	gshiftleft(Q, N);
	itog(1, t1);
	addg(t1, N);
	pr[0] = 2;
	for (k=0, npr=1;; k++)
	{
		if (primeq(3+2*k))
		{
			pr[npr++] = 3+2*k;
			if (npr >= NUM_PRIMES)
				break;
		}
	}

	set_random_seed();
	pass = 0;
	while (++pass)
	{
        if(test_mode) {
        	switch(atoi(argv[1])) {
/* Starting seed will factor F7 in stage 2. */
           		case 7:  seed = 442057745 + pass -1; B = 100000L; stage_flag = 2; break;
           		case 8:  seed = 1858194296 + pass -1; B = 100000L; stage_flag = 2; break;
           		case 9:  seed = 361900231 + pass -1; B = 100000L; stage_flag = 1; break;
           		case 10:  seed = 1587102819 + pass -1; B = 100000L; stage_flag = 1; break;
           		case 11:  seed = 2034710845 + pass -1; B = 100000L; stage_flag = 1; break;
				case 12: seed = 1266569162 + pass - 1; stage_flag = 2; break;
/* Starting seed will factor F13 in stage 1. */
           		case 13: seed = 8020345 + pass-1; B = 70000L; stage_flag = 1; break; 
           		case 15: seed = 125546653 + pass-1; B = 120000L; stage_flag = 2; break; 
           		case 16: seed = 253301772 + pass-1; B = 600000L; stage_flag = 2; break; 
           		case 17: seed = 453575255 + pass-1; B = 100000L; stage_flag = 2; break; 
           		case 18: seed = 1608584961 + pass-1; B = 100000L; stage_flag = 2; break; 
          		default: seed = psi_rand(); B = 100000L; stage_flag = 0; break;
	    	} 
            C = 50*B;
		    if(stage_flag != 0) printf("TESTING: This test should find a factor in stage %d.\n", stage_flag);
		}  else seed = psi_rand();
        choose12(xr, zr, seed, An, Ad, N);
        printf("Selecting elliptic curve %d seed s = %d; B = %d, C = %d:\n", pass, seed, B, C);
		fflush(stdout);
		cnt = 0;
		nshorts = 1;
		count = 0;
		for (j=0;j<nshorts;j++)
		{
			ell_mul(xr, zr, (1 << 16), An, Ad, N);
			ell_mul(xr, zr, 3*3*3*3*3*3*3*3*3*3*3, An, Ad, N);
			ell_mul(xr, zr, 5*5*5*5*5*5*5, An, Ad, N);
			ell_mul(xr, zr, 7*7*7*7*7*7, An, Ad, N);
			ell_mul(xr, zr, 11*11*11*11, An, Ad, N);
			ell_mul(xr, zr, 13*13*13*13, An, Ad, N);
			ell_mul(xr, zr, 17*17*17, An, Ad, N);
		}
		k = 19;
		while ((unsigned int)k<B)
		{
			if (isprime(k))
			{
				ell_mul(xr, zr, k, An, Ad, N);
				if (k<500)
					ell_mul(xr, zr, k, An, Ad, N);
				if (cnt++ %1000==0)
					dot();
			}
			k += 2;
		}
		count = 0;

		gtog(zr, gg);
		gcdg(N, gg);
		if((!isone(gg))&&(bitlen(gg)>limitbits))
		{
			fprintf(stdout,"\n!{s, B, C} = {%d, %u, %u};  ", seed, B, C);
			gwriteln(gg, stdout);
			fflush(stdout);
            if(test_mode) if(stage_flag == 1) exit(0);
		}
		else
		{
			printf("\n");
			fflush(stdout);
		}
        
        if(test_mode) if(stage_flag == 1) exit(0);

		printf("\nCreating second stage table...\n");
		fflush(stdout);
		k = (int)(B/D);
		gtog(xr, xb[0]); gtog(zr, zb[0]);
		ell_mul(xb[0], zb[0], k*D+1 , An, Ad, N);
		gtog(xr, xb[D+1]); gtog(zr, zb[D+1]);
		ell_mul(xb[D+1], zb[D+1], (k+2)*D+1 , An, Ad, N);

		for (j=1; j <= D; j++)
		{
			gtog(xr, xb[j]);
			gtog(zr, zb[j]);
			ell_mul(xb[j], zb[j], 2*j , An, Ad, N);
			gtog(zb[j], xzb[j]);
			special_mulg(xb[j], xzb[j]);
			s_modg(N, xzb[j]);
			if(j%10 == 0) dot();
		}
		printf("\nCommencing second stage...\n");
		fflush(stdout);
		count = 0;
		itog(1, gg);

		while (1)
		{
			gtog(zb[0], xzb[0]);
			special_mulg(xb[0], xzb[0]);
			s_modg(N, xzb[0]);
			special_mulg(zb[0], gg);
			s_modg(N,gg); /* Accumulate. */

/* Next, process sieve. */
/* Next, initialize sieving offsets. */
        	parr[0] = 0;
        	for(j=1; j < D; j++) parr[j] = 1;
            m = (k+2)*D;
			m1 = (k*D + 3);
			for(i = 1; ; i++) {
                pri = pr[i];
				if(pri*pri > m) break;
				m2 = -(m1 % pri);
			   if(m2 < 0)  m2 += pri;
			   pjump[i] = m2;
		    }
			iend = i;
            m2 = 2*D-2;
		    for(i = 1; i < iend; i++) {
				m1 = pr[i];
				for(m = pjump[i]+1; m < m2; m += m1) {
						if(m & 1) parr[(m+1)>>1] = 0;
				}
			}
			for (j = 1; j < D; j++)
			{
/* Integrity test for sieve: */
				if(!parr[j])
					continue;
				/* Next, accumulate (xa - xb)(za + zb) - xa za + xb zb. */
				gtog(xb[0], t1);
				subg(xb[j], t1);
				gtog(zb[0], t2);
				addg(zb[j], t2);
				special_mulg(t1, t2);
				s_modg(N, t2);
				subg(xzb[0], t2);
				addg(xzb[j], t2);
				special_mulg(t2, gg);
				s_modg(N, gg);
				if ((++count)%1000==0)
					dot();
			}

			k += 2;
			if (((unsigned long)k)*D >= C)
				break;
			gtog(xb[D+1], xs);
			gtog(zb[D+1], zs);
			ell_odd(xb[D], zb[D], xb[D+1], zb[D+1], xb[0], zb[0], N);
			gtog(xs, xb[0]);
			gtog(zs, zb[0]);
		}

		gcdg(N, gg);
		if ((!isone(gg))&&(bitlen(gg)>limitbits))
		{
			fprintf(stdout,"\n!!{s, B, C} = {%d, %u, %u};  ", seed, B, C);
			gwriteln(gg, stdout);
			fflush(stdout);
            if(test_mode) if(stage_flag == 2) exit(0);
			continue;
		}
        if(test_mode) if(stage_flag == 2) exit(0);
		printf("\n"); fflush(stdout);
	}
	return 0;
}
