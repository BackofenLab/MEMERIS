/*
 * $Id: mast.h,v 1.2 2005/10/20 00:20:52 tbailey Exp $
 * 
 * $Log: mast.h,v $
 * Revision 1.2  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2005/07/29 18:43:40  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef MAST_H 
#define MAST_H 

#include "meme.h"
#include "macros.h"
#include "general.h"
#include "banner.h"
#include "pssm-distr.h"
#include <sys/types.h>
#include <sys/stat.h>

/*
	Macro definitions
*/

/* compute the p-value of value of n samples of an extreme value whose 
   p-value is known and return it in r
*/ 
#define EV(pvalue, n, r) \
        (r) = (pvalue) * (n);		/* good approximation if small */ \
        if (r > .001) (r) = 1.0 - pow((1.0 - (pvalue)), (double)(n));

/*
	Data definitions 
*/

/* scores for each position in a sequence */
typedef struct score {
  double score;				/* the score */
  BOOLEAN ic;				/* score of inverse complement strand */
} SCORE;

/*
	Procedures
*/

extern SCORE **score_sequence(
  STYPE stype,		/* handling of different DNA strands */
  BOOLEAN xlate_dna,	/* database is DNA and motifs protein */
  char *sequence,	/* the sequence to score */
  long length,		/* length of the sequence */
  int nmotifs, 		/* number of motifs */
  LO *los[]   		/* array of pointers to log-odds matrices */
);

extern TILING score_tile_diagram(
  char *sequence,			/* sequence to score and tile */
  long length,				/* length of sequence */
  LO *los[],				/* array of pointers to lo matrices */
  int nmotifs,				/* number motifs read */
  BOOLEAN dna,                          /* database is DNA */
  STYPE stype,				/* handling of different strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN best_motifs,                  /* show only best motifs in diagrams */
  BOOLEAN print_p,			/* print p-value in block if TRUE */
  double **pv,				/* p-value tables for each motif */
  double m_thresh,			/* maximum motif p-value to print */
  double w_thresh,			/* maximum motif p-value-- weak hits */
  BOOLEAN use_seq_p, 			/* use sequence not position p-values */
  BOOLEAN hit_list                      /* create hit list instead of diagram */
);

extern void print_diagram(
  char *dia,				/* motif diagram string */
  char *hdr,				/* prefix for each line of diagram */
  FILE *file				/* destination file */
);

extern void make_block(
  int m,			/* motif number */
  char *strand,			/* strand */
  int f,			/* frame number; f=0 not translating DNA */
  double thresh,		/* strong motif threshold */
  double p,			/* p-value */
  BOOLEAN print_p,		/* print p-value in block if TRUE */
  char *block			/* put block string here */
);	

extern double qfast(
  int n,			/* number of random variables in product */
  double k			/* product of random variables */
);

extern double *init_mast_background(
  char *bfile,				/* name of background file */
  char *alphabet,			/* motif alphabet */
  STYPE stype,				/* handling of DNA strands */
  BOOLEAN translate_dna,		/* DNA sequences and protein motifs */
  char **blast_alphabet,		/* corresponding BLAST alphabet */
  int *p[MAXASCII]			/* permutation/substitution matrix */
);

extern void free_tiling(
  TILING tiling
); 

extern char *create_diagram(
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* treatment of strands of DNA */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN best_motifs, 			/* diagrams have only best motif */
  BOOLEAN print_p,			/* print p-value in block */
  double thresh,			/* strong hit threshold */
  int nmotifs,				/* number of motifs */
  LO *los[],				/* array of pointers to lo matrices */
  long length,				/* length of sample */
  BOOLEAN hit_list,			/* create hit list instead of diagram */
  TILING tiling				/* tiling of sequence */
);

double *get_seq_comp(
  BOOLEAN xlate_dna,			/* translate DNA */
  char *sequence,			/* ASCII sequence */
  int alen				/* length of alphabet */
);

#endif

