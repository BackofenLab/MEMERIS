/*
 * $Id: star.h,v 1.1.1.1 2005/07/29 19:11:34 nadya Exp $
 * 
 * $Log: star.h,v $
 * Revision 1.1.1.1  2005/07/29 19:11:34  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef STAR_H
#define STAR_H

/* range of scores in pairwise score matrix */
#define PAIR_RANGE 100

extern void init_star_distr(
  THETA logodds,                        /* pairwise score table (logodds) */
  int N,				/* maximum sites in star alignments */
  int alen,				/* length of alphabet */
  double *back				/* letter distribution; include p(X)=0*/
);

extern THETA init_pairwise_matrix(
  MAP_TYPE map_type,			/* type of mapping */
  double map_scale,			/* scale of mapping */
  int alength, 				/* length of alphabet */
  double *back				/* background frequencies */
);

extern double get_star_logpv(
  double score,				/* star alignment score */
  double wN				/* weighted number of sites in algnmt.*/
);

#endif
