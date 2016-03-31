/*
 * $Id: gamma.h,v 1.1.1.1 2005/07/29 18:38:01 nadya Exp $
 * 
 * $Log: gamma.h,v $
 * Revision 1.1.1.1  2005/07/29 18:38:01  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/* 
	To compute the gamma function and the incomplete gamma function. 
*/

#ifndef GAMMA_H
#define GAMMA_H

#define ITMAX 100
#define EPS 1e-20

/* the incomplete gamma Q function */
#define GAMMAQ(a,x)	( ((x) < (a)+1) ? (1.0 - gser((a),(x))) : gcf((a),(x)))

extern double gser(
  double a, 
  double x
); 

extern double gcf(
  double a, 
  double x
);

extern double loggamma(
  double x
);
 
#endif
