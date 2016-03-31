/*
 * $Id: background.h,v 1.1.1.1 2005/07/29 18:34:37 nadya Exp $
 * 
 * $Log: background.h,v $
 * Revision 1.1.1.1  2005/07/29 18:34:37  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/**********************************************************************/
/*
	MEME background models
*/
/**********************************************************************/

#ifndef BACKGROUND_H
#define BACKGROUND_H

/* maximum size of background model */
#define MAX_BACK_SIZE  19530

/* compute the log probability of a substring of a sequence given the log
   cumulative probability in lcb:
   log Pr(S_{i,...,i+w-1} | H_0) 
*/
#define Log_back(lcb, i, w) (lcb[(i)+w] - lcb[i])

extern double *read_markov_model(
  char *pfile,					/* name of probability file */
  double *freq,					/* letter frequencies */
  char *alpha,					/* alphabet expected */
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  BOOLEAN rc,					/* average reverse complements*/
  int *order					/* order of model read */
);

extern double log_cum_back(
  char *seq,					/* the sequence */
  double *a_cp,					/* the conditional probs */
  int order,					/* the order of the model */
  double *logcumback				/* the result */
);

#endif

