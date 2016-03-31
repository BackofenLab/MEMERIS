/*
 * $Id: dpalign.h,v 1.1.1.1 2005/07/29 18:37:27 nadya Exp $
 * 
 * $Log: dpalign.h,v $
 * Revision 1.1.1.1  2005/07/29 18:37:27  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef DPALIGN_H
#  define DPALIGN_H

extern char *dp_align(
  int alength,					/* length of alphabet */
  int n,					/* length of logodds matrix */
  LOGODDS logodds,				/* logodds matrix */
  char *seq1,					/* consensus of logodds */
  int m,					/* length of sequence */
  char *eseq,					/* integer encoded sequence */
  double wg,					/* gap cost (initialization) */ 
  double ws,					/* space cost (extension) */ 
  BOOLEAN endgaps 				/* penalize end gaps if TRUE */
);

extern char **dp_multi_align(
  P_PROB sites,                                 /* the sites */
  int nsites,                                   /* the number of sites */
  int w,                                        /* width of sites */
  int flank,					/* add flank cols on lft+rgt */
  DATASET *dataset                              /* the dataset */
);

extern int g_align(
  char **ma,				/* the multiple alignment */
  int nseq,				/* number of sequences in alignment */
  int ncol,				/* number of columns in alignment */
  int left,				/* leftmost allowed end */
  int right,				/* rightmost allowed end */
  int minw,				/* minimum width */
  int *off,				/* offset of g-alignment rel. to left */
  int *w				/* width of g-alignment */
);

#endif
