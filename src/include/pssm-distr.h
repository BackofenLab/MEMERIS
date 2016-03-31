/*
 * $Id: pssm-distr.h,v 1.2 2005/10/25 21:23:44 nadya Exp $
 * 
 * $Log: pssm-distr.h,v $
 * Revision 1.2  2005/10/25 21:23:44  nadya
 * add Id and Log lines for keeping track of changes
 * change C++-style comments to C-style
 *
 *
 */


#ifndef pssm_distr_h
#define pssm_distr_h

extern double *calc_pssm_cdf(
  int w,                /* width of PSSM */
  int alen,             /* length of alphabet */
  int range,            /* largest value in PSSM */
  double **pssm,        /* scaled, integer PSSM: pssm[i][j] is score for 
                           j_th letter in i_th column of motif 
                           entries in PSSM are in range [0..R */
  double *prob          /* 0-order Markov background mode */
);

#endif
