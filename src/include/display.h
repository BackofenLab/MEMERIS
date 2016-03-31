/*
 * $Id: display.h,v 1.1.1.1 2005/07/29 18:37:03 nadya Exp $
 * 
 * $Log: display.h,v $
 * Revision 1.1.1.1  2005/07/29 18:37:03  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */


#ifndef DISPLAY_H
#  define DISPLAY_H

extern void print_results(
  DATASET *dataset,		/* the dataset */
  DATASET *neg_dataset,		/* negative examples */
  MODEL *model,			/* the model */
  MODEL *neg_model, 		/* negative model */
  CANDIDATE *candidates  	/* candidate models found */
);

extern void print_theta(
  int imotif,		/* motif number */
  int format,		/* 1 = floating point
              		   2 = integer */ 
  int nsites,           /* number of sites (discrete) */
  THETA theta,		/* theta */
  int w,		/* width of motif */
  double log_ev,	/* log motif E-value */
  char *str_space,	/* space for printing strand direction */
  DATASET *dataset,	/* the dataset */ 
  FILE *outfile	 	/* NULL stdout; otherwise, print to split file */
);

extern void print_zij(
  DATASET *dataset,			/* the dataset */
  MODEL *model				/* the model */
);

extern void print_wij(
  DATASET *dataset		/* the dataset */
);

extern char *get_consensus(
  THETA theta,			/* motif theta */
  int w,			/* width of motif */
  DATASET *dataset,		/* the dataset */
  int N,			/* number of letters for each position */
  double min_prob		/* minimum cumulative prob for N letters */
);

extern void print_command_summary(
  MODEL *model,                 /* the model */
  DATASET *dataset              /* the dataset */
);

extern void print_dataset_summary (
  DATASET *dataset 		/* the dataset */
);

extern void print_summary(
  MODEL *model,				/* the model */
  DATASET *dataset,                     /* the dataset */
  LO **los,                             /* the LO structures */
  int nmotifs,                          /* number of motifs */ 
  double **pv,                          /* p-value of score distribution */
  FILE *outfile                         /* where to print */
);

extern void print_meme_doc(
  char *meme_directory                  /* meme source directory */
);

#endif

