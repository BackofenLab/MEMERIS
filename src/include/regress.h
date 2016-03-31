/*
 * $Id: regress.h,v 1.1.1.1 2005/07/29 19:09:59 nadya Exp $
 * 
 * $Log: regress.h,v $
 * Revision 1.1.1.1  2005/07/29 19:09:59  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1995, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/
 
#ifndef REGRESS_H 
#define REGRESS_H

extern double regress(
  int n,                        /* number of points */
  double *x,                    /* x values */
  double *y,                    /* y values */
  double *m,                    /* slope */
  double *b                     /* y intercept */
);

#endif
 
