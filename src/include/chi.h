/*
 * $Id: chi.h,v 1.1.1.1 2005/07/29 18:35:52 nadya Exp $
 * 
 * $Log: chi.h,v $
 * Revision 1.1.1.1  2005/07/29 18:35:52  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef CHI_H
#define CHI_H

extern double chiq(
  double chisq,			/* chi-square variable */
  double v 			/* degrees of freedom */
);

extern double chi(
  double nu,    /* degrees of freedom */
  double alpha
);

#endif
