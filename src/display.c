/*
 * $Id: display.c,v 1.3.4.1 2006/01/24 20:44:08 nadya Exp $
 * 
 * $Log: display.c,v $
 * Revision 1.3.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.3  2005/10/25 19:00:37  nadya
 * change c++ style comment to proper c
 *
 * Revision 1.2  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.1.1.1  2005/07/29 23:35:53  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME							       *
*	Copyright 1994-2006, The Regents of the University of California *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/
/**********************************************************************/
/*
	Display routines for the results of MEME
*/
/**********************************************************************/
/* 7-23-99 tlb; replace nsites() calls with model->nsites */
/* 7-16-99 tlb; move SITE to meme.h, get_sites to meme_util.c */
/* 7-14-99 tlb; change simplified prob. matrix to be observed frequencies,
  consensus to be based on observed frequencies */
/* 4-7-99 tlb; fix bug in get_sites: setting of best site in oops and zoops */

#include <meme.h>
#include <mast.h>

static char *yesno[] = {"no", "yes"};
static char *stars = NULL;

/* number of different integral score values for computing p-values */
#define RANGE 100

/* stepsize and size of smoothing window for get_q */
#define NSTEPS 100
/*#define WINDOW (NSTEPS/10)*/
#define WINDOW 0

/* round double to integer; round up if half way between */
/* only works for positive numbers */
#ifndef NINT
  #define NINT(x) ((int) ((x)+0.5))
#endif

/* maximum number of levels in consensus sequence */
#define MAXDEPTH ((int) (1.0/MINCONS))

/* encode and decode a (seq #, offset) pair as one integer */
/* this will bomb if seq. length > MAX_SEQ_LEN or
   if number of sequences * MAX_SEQ_LEN exceeds size of int 
*/
#define MAX_SEQ_LEN 10000
#define ENCODE(I,J) (I) * MAX_SEQ_LEN + (J);
#define DECODE(N, SEQ, OFF) {\
  (OFF) = (N) % MAX_SEQ_LEN;\
  (SEQ) = ((N) - (OFF)) / MAX_SEQ_LEN;\
}

/* distance to indent start of IC histogram, consensus and simplified motif */
#define IND 13
#define IND2 6

/* record describing accuracy of motif */
typedef struct {
  double thresh;	/* optimal threshold */
  int err;		/* classification errors using thresh */
  double roc;		/* ROC */
} ACCURACY;

/* sortable sequence score record */
typedef struct {
  double score;		/* sequence score */
  int class;		/* class of sequence,  1=pos, 0=neg */
  char *id;
} SORTED_SCORE;

/* local functions */
static void print_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  int format,				/* 0=BLOCKS; 1=FASTA */
  char *com,				/* comment to append */
  FILE *outfile				/* where to print */
);
static void print_log_odds(
  int imotif,			/* motif number */ 
  DATASET *dataset,		/* the dataset */
  int w,			/* width of motif */
  BOOLEAN pair,			/* double matrix if true */
  LOGODDS logodds,		/* log-odds matrix */
  double bayes,  		/* threshold */
  double log_ev			/* log E-value of motif */
);
static void print_entropy(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  char *str_space,		/* space for printing strand direction */
  FILE *outfile			/* stream for output */
);
static void print_candidates(
  CANDIDATE *candidates, 	/* list of candidate models */
  DATASET *dataset,		/* the dataset */
  int max_w			/* maximum width for motifs */
);
static ACCURACY *get_thresh(
  int w,				/* width of motif */
  LOGODDS logodds1,			/* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2,			/* log-odds matrix: LOG2(m_ij/a_ij) */
  DATASET *pos,				/* positive examples */
  DATASET *neg,				/* negative examples */
  BOOLEAN print_scores			/* print sorted scores */
);
static double meme_score_sequence(
  char *eseq,		/* integer-coded sequence to score */
  int length,		/* length of the sequence */
  int w, 		/* width of motif */
  LOGODDS logodds1, 	/* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2  	/* log-odds matrix: LOG2(m_ij/n_ij) */
);
static int s_compare(
  const void *v1,
  const void *v2
);
static double get_q(
  int nsteps,					/* try nsteps from 0 to 1 */
  int window,					/* smoothing window radius */
  int w,					/* width of motif */
  THETA theta,					/* motif theta */
  THETA neg_theta, 				/* anti-motif theta */
  double *back,					/* background motif */
  DATASET *dataset,				/* the dataset */
  DATASET *neg_dataset,				/* negative examples */
  char *str_space		/* space for printing strand direction */
);
static LO *create_lo(
  int name,				/* name of motif */
  int w, 				/* width of motif */
  int alen,				/* length of alphabet */
  LOGODDS logodds,			/* single-letter logodds matrix */
  BOOLEAN pair,
  double threshold			/* Bayes optimal threshold */
);
static void score_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  LO *lo,				/* LO structure */
  double *pv 				/* p-values for scores of this motif */
);
static void print_site_diagrams(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  int nmotifs,				/* number of motifs in los */
  LO *los[MAXG],			/* logodds structure for motif */
  FILE *outfile				/* where to print */
);
static void align_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  LO *lo,				/* LO structure */
  double *pv,				/* pvalues for scores of this motif */
  FILE *outfile				/* stream for output */
);
static void print_block_diagrams(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  LO *los[MAXG],			/* logodds structures for motifs */
  int nmotifs, 				/* number of motifs */
  double *pv[MAXG],			/* p-value tables for each motif */
  FILE *outfile
);

/**********************************************************************/
/*
	Print the results of EM
*/
/**********************************************************************/
extern void print_results(
  DATASET *dataset,		/* the dataset */
  DATASET *neg_dataset,		/* negative examples */
  MODEL *model,			/* the model */
  MODEL *neg_model, 		/* negative model */
  CANDIDATE *candidates 	/* candidate models found */
)
{
  int i;
  int max_w = model->max_w;			/* maximum width tried */
  int nstrands = model->invcomp ? 2 : 1;	/* # of strands to use */
  int w = model->w;				/* width of last component */
  int nsites_dis = model->nsites_dis;	 	/* # of sites after discretiz.*/
  double m1, e1, m2, e2;			/* for printing significance */
  char *str_space = (nstrands == 1) ? "" : "       ";
  BOOLEAN pair = neg_dataset ? neg_dataset->negtype == Pair : FALSE;
  BOOLEAN blend = neg_dataset ? neg_dataset->negtype == Blend : FALSE;
  THETA theta = model->theta;
  THETA obs = model->obs;
  THETA neg_theta = neg_model ? neg_model->theta : NULL;
  THETA neg_obs = neg_model ? neg_model->obs : NULL;
  double lambda = model->lambda;
  LOGODDS logodds, logodds1;
  ACCURACY *acc=NULL;
  char *cons = model->cons;
  int alength = dataset->alength;
  double *back = dataset->back;
  int imotif = model->imotif-1;			/* index of motif */
  double thresh;				/* Bayes optimal threshold */
  FILE *outfile = stdout;

  /* get entropy and E-value of motif */
  calc_entropy(model, dataset);

  /* get p-value and E-value of motif */
  exp10_logx(model->logpv/log(10.0), m1, e1, 1);
  exp10_logx(model->logev/log(10.0), m2, e2, 1);

  /* print the significant models */
  if (VERBOSE) {
    print_candidates(candidates, dataset, max_w);
  }

  /* print the results for the model as a whole */
  if (pair || blend) {			/* negative motifs */
    int n_nsites_dis = neg_model->nsites_dis;
    printf( "\n\n%s\nMOTIF%3d\twidth = %3d\tsites = %4d\tnegative sites = %4d",
      stars, model->imotif, w, nsites_dis, n_nsites_dis);
  } else {
    printf("\n\n%s\nMOTIF%3d\twidth = %4d   sites = %3d   ",
      stars, model->imotif, w, nsites_dis);
  }
  printf("llr = %.0f   E-value = %3.1fe%+04.0f\n%s\n", model->llr, m2, e2, stars);

  if (VERBOSE) {
    printf("p-value = %3.1fe%+04.0f   E-value = %3.1fe%+04.0f\n%s\n", 
      m1, e1, m2, e2, stars);
  }

  /* print results for motif */
  if (VERBOSE) {
    char *cons0 = model->cons0;
    printf("\n(best) %s --> %s\n", cons0, cons);
    printf("(best) w %3d nsites %5.1f lambda %8.7f IC/col %6.3f\n",
      w, lambda*wps(dataset, w), lambda, model->rel);
    printf("\n(best) IC %6.3f\n\n", w * model->rel);
  }

  /* 
    print motif description
  */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d Description\n", model->imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);

  /* print the simplified (1-digit) motif */
  /* make the log-odds matrices */
  /* calculate the optimal threshold (min classification error or Bayes' */
  /* print the IC of the motif as a bar graph */
  if (blend) {			/* single blended motif */
    /* get mixing parameter */
    double q = get_q(NSTEPS, WINDOW, w, theta, neg_theta, back,
      dataset, neg_dataset, str_space);
    /* get threshold and performance of combined motif */
    logodds = make_log_odds(theta, neg_theta, back, q, w, alength);
    acc = get_thresh(w, logodds, NULL, dataset, neg_dataset, FALSE);
    thresh = acc->thresh;
    myfree(acc);
    printf("\t\tPOSITIVE MOTIF\n");
    print_theta(0, 2, model->nsites_dis, obs, w, model->logev, 
      str_space, dataset, stdout);
    printf("\t\tNEGATIVE MOTIF\n");
    print_theta(0, 2, neg_model->nsites_dis, neg_obs, w, model->logev, 
      str_space, dataset, stdout);
    print_entropy(neg_model, dataset, str_space, stdout);
  } else if (pair) {			/* pair of motifs, pos and neg */
    logodds = make_log_odds(theta, NULL, back, 0, w, alength);
    logodds1 = make_log_odds(neg_theta, NULL, back, 0, w, alength);
    thresh = LOG2((1-lambda)/lambda);	/* Bayes' threshold */
    printf("\t\tPOSITIVE MOTIF\n");
    print_theta(0, 2, model->nsites_dis, obs, w, model->logev, 
      str_space, dataset, stdout);
    print_entropy(model, dataset, str_space, stdout);
    printf("\t\tNEGATIVE MOTIF\n");
    print_theta(0, 2, neg_model->nsites_dis, neg_obs, w, model->logev, 
      str_space, dataset, stdout);
    print_entropy(neg_model, dataset, str_space, stdout);
    myfree(logodds1);
    /* combine the matrices by appending P/N to P/B */
    logodds1 = make_log_odds(theta, neg_theta, NULL, 1, w, alength);
    Resize(logodds, 2*w, double *);
    for (i=0; i<w; i++) logodds[w+i] = logodds1[i];
    myfree(logodds1);
  } else {				/* no negative motifs */
    logodds = make_log_odds(theta, NULL, back, 0, w, alength);
    thresh = LOG2((1-lambda)/lambda);	/* Bayes' threshold */
    print_theta(0, 2, model->nsites_dis, obs, w, model->logev, 
      str_space, dataset, stdout);
    print_entropy(model, dataset, str_space, stdout);
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n"); 

  /* create an LO structure and store it in the local array */
  los[imotif] = create_lo(model->imotif, w, alength, logodds, pair, thresh);
  
  /* create a table of p-values and store it in the array */
  los[imotif]->alen--;			/* ignore X since freq isn't set right*/
  /*pv[imotif] = calc_cdf(los[imotif], RANGE, dataset->back);*/
  pv[imotif] = calc_pssm_cdf(los[imotif]->w, los[imotif]->alen, RANGE, los[imotif]->logodds, dataset->back);

  los[imotif]->alen++;

  /* score the sites and sort by position p-value */
  score_sites(dataset, model, los[imotif], pv[imotif]);

  /* print alignment of the sites */
  align_sites(dataset, model, los[imotif], pv[imotif], stdout);

  /* print diagrams of the sites */
  print_site_diagrams(dataset, model, model->imotif, los, stdout);

  /* print the sites "making up" the model */
  if (pair || blend) {
    print_sites(dataset, model, PRINT_FASTA, " (positive motif)", stdout);
    print_sites(neg_dataset, neg_model, PRINT_FASTA, " (negative motif)", 
      stdout);
  } else {
    print_sites(dataset, model, PRINT_FASTA, "", stdout);
  }

  /* print the logodds matrix */
  print_log_odds(model->imotif, dataset, w, pair, logodds, thresh,model->logev);

  /* print the observed frequency matrix */
  print_theta(model->imotif, 1, model->nsites_dis, obs, w, model->logev, 
    str_space, dataset, stdout);
  if (pair || blend) 
    print_theta(model->imotif, 1, neg_model->nsites_dis, neg_obs, w, 
      model->logev, str_space, dataset, stdout);

  /* display elapsed time */
  printf("\n\n\nTime %5.2f secs.\n\n", myclock()/1E6);

  /* print line of stars */
  fprintf(stdout, "%s\n", stars);

  /* flush */
  fflush(stdout);

} /* print_results */

/**********************************************************************/
/*
	create_lo

	Create an an LO structure from the logodds matrix;
	include 'X' in the alphabet.
*/
/**********************************************************************/
static LO *create_lo(
  int name,				/* name of motif */
  int w, 				/* width of motif */
  int alen,				/* length of alphabet */
  LOGODDS logodds,			/* single-letter logodds matrix */
  BOOLEAN pair,
  double threshold			/* Bayes optimal threshold */
)
{
  int i, j;
  LO *lo = NULL;			/* the LO structure */
  
  /* create a logodds structure */
  Resize(lo, 1, LO);

  /* initialize it */
  lo->alen = alen+1;			/* add 'X' column */
  lo->w = lo->ws = w;
  lo->pair = pair;
  lo->name = name;
  lo->thresh = threshold;

  /* make a copy of the logodds matrix and scale it to [0..range] */
  create_2array(lo->logodds, LOGODDSB, w, lo->alen);
  for (i=0; i<w; i++) for(j=0; j<lo->alen; j++) lo->logodds(i,j) = logodds(i,j);
  scale_lo(&lo, 1, RANGE);				/* scale */
  make_double_lo(&lo, 1);	/* make a double-letter logodds matrix */

  return(lo);
} /* create_lo */

/**********************************************************************/
/*
	print_block_diagrams

	Tile the dataset sequences with the motifs in los[] and print
	the block diagrams with the p-value of the product of p-values.
*/
/**********************************************************************/
static void print_block_diagrams(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  LO *los[MAXG],			/* logodds structures for motifs */
  int nmotifs, 				/* number of motifs */
  double *pv[MAXG],			/* p-value tables for each motif */
  FILE *outfile
)
{
  int i;
  BOOLEAN dna = dataset->dna;		/* dataset is DNA */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand */
  STYPE stype = dna ? (invcomp ? Combine : Norc) : Protein;
  BOOLEAN xlate_dna = FALSE;		/* don't translate DNA */
  BOOLEAN best_motifs = FALSE;		/* use all motifs */
  double m_thresh = 1e-4;		/* show motifs over p-value 0.0001 */
  double w_thresh = m_thresh;		/* show strong motifs only */
  BOOLEAN use_seq_p = FALSE;		/* use postion p-values */
  char *f = "%-*s%s %8s  %s\n";		/* format */

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tCombined block diagrams:");
  fprintf(outfile, " non-overlapping sites with p-value < %6.4f\n", m_thresh);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "COMBINED P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  for (i=0; i<dataset->n_samples; i++) {
    SAMPLE *s = dataset->samples[i];
    char *name = s->sample_name;
    int lseq = s->length;
    char *sequence = s->seq;
    TILING tiling = score_tile_diagram(sequence, lseq, los, nmotifs,
      dna, stype, xlate_dna, best_motifs, TRUE, pv, m_thresh, w_thresh, 
      use_seq_p, FALSE);
    fprintf(outfile, "%-*.*s %16.2e  %s\n", MSN, MSN, name, tiling.pv,
      tiling.diagram);
    free_tiling(tiling);
  } /* sequence */

  /* print a final line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");
  
} /* print_block_diagrams */

/**********************************************************************/
/*
	print_log_odds

	Print the log-odds matrix 
*/
/**********************************************************************/
static void print_log_odds(
  int imotif,			/* motif number */ 
  DATASET *dataset,		/* the dataset */
  int w,			/* width of motif */
  BOOLEAN pair,			/* double matrix if true */
  LOGODDS logodds,		/* log-odds matrix */
  double bayes,  		/* threshold */
  double log_ev			/* log E-value of motif */
)
{
  int i, j;
  int alength = dataset->alength;	/* length of alphabet */
  int n = wps(dataset, w);		/* weighted possible sites */
  char *type = pair ? "pair" : "";	/* type of matrix */
  double m1, e1;			/* for printing significance */
  FILE *outfile = stdout;

  /* get E-value of motif */
  exp10_logx(log_ev/log(10.0), m1, e1, 1);

  /* double w if printing a matrix pair */
  if (pair) w *=2;

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\tMotif %d position-specific scoring matrix\n", imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n");
  fprintf(outfile, 
    "log-odds matrix: alength= %d w= %d n= %d bayes= %g E= %3.1fe%+04.0f %s\n", 
    alength, w, n, bayes, m1, e1, type);

  for (i=0; i < w; i++) {		/* site position */
    for (j=0; j < alength; j++) {	/* letter */
      fprintf(outfile, "%6d ", NINT(100*logodds(i,j)));
    }
    fprintf(outfile, "\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");
} /* print_log_odds */

/**********************************************************************/
/*
	print_entropy	

	Display the relative entropy of each column of the motif
	as a bar graph.
*/
/**********************************************************************/
static void print_entropy(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  char *str_space,		/* space for printing strand direction */
  FILE *outfile			/* stream for output */
)
{
  int i, j;
  int w = model->w;				/* width of motif */
  THETA obs = model->obs;			/* observed frequencies */
  double *rentropy = model->rentropy;		/* IC of each column */
  double ic = model->rel * w;			/* motif information content */
  double *back = dataset->back;			/* background model */
  int alength = dataset->alength;
  char icstring[15];				/* print string for ic */
  char *consensus;				/* consensus strings */
  double min_freq;				/* minimum background freq */
  double maxre;					/* maximum relative entropy */
  int nsteps;					/* number of steps histogram */

  /* get minimum background frequency and maximum relative entropy */
  for (i=0, min_freq=1; i<alength; i++)  
    if (back[i] < min_freq) min_freq = back[i];
  maxre = -LOG2(min_freq);		/* maximum relative entropy */
  
  /* create string containing IC */
  sprintf(icstring, "(%.1f bits)", ic); 

  /* print the relative entropy of each column as a bar graph */
  nsteps = 10;
  for (i=0; i<nsteps; i++) {
    double level = maxre - (i * maxre/nsteps);
    fprintf(outfile, (i==0 ? "%*.*s %*.1f " : "%-*.*s %*.1f "), IND, IND, 
      (i==0 ? "bits" : i==4 ? "Information" : i==5 ? "content" : 
        i==6 ? icstring : ""), IND2, level);
    for (j=0; j<w; j++) {
      if (NINT(nsteps * rentropy[j] / maxre) >= nsteps-i) {
        fputc('*', outfile);
      } else {
        fputc(' ', outfile);
      }
    }
    fputc('\n', outfile);
  }
  fprintf(outfile, "%-*.*s %*.1f ", IND, IND, "", IND2,0.0);
  for (i=0; i<w; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

  /* get and print the consensus sequences */
  consensus = get_consensus(obs, w, dataset, MAXDEPTH, MINCONS);
  for (i=0; i<MAXDEPTH && i<alength; i++) {/* print next levels of consensus */
    fprintf(outfile, "%-*.*s %*.0s %*.*s\n", IND, IND, 
      (i==0 ? "Multilevel" : i == 1 ? "consensus" : i == 2 ? "sequence" : ""), 
      IND2, "", w, w, consensus+(i*w));
  }

  /* free up space */
  myfree(consensus);
} /* print_entropy */

/**********************************************************************/
/*
	print_theta

		format=1		floating point; pos x letter
		format=2		1 digit; letter x pos 

	Print the probability array.
*/
/**********************************************************************/
extern void print_theta(
  int imotif,		/* number of motif */
  int format,		/* 1 = floating point
              		   2 = integer */ 
  int nsites,		/* number of sites (discrete) */
  THETA theta,		/* theta */
  int w,		/* width of motif */
  double log_ev,	/* log motif E-value */
  char *str_space,	/* space for printing strand direction */
  DATASET *dataset,	/* the dataset */ 
  FILE *outfile	 	/* file to print to */
)
{
  int i, j;
  int alength = dataset->alength;

  if (format == 1) {
    double e1, m1;
    exp10_logx(log_ev/log(10.0), m1, e1, 3);
    for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
    fprintf(outfile, "\n\tMotif %d position-specific probability matrix\n", 
      imotif);
    for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
    fprintf(outfile,
      "\nletter-probability matrix: alength= %d w= %d nsites= %d E= %3.1fe%+04.0f ",
      alength, w, nsites, m1, e1);
    fprintf(outfile, "\n");
    for (i=0; i < w; i++) {
      for (j=0; j < alength; j++) {
	fprintf(outfile, "%9.6f ", theta(i, j));
      }
      fprintf(outfile, "\n");
    }
    for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
    fprintf(outfile, "\n\n");

  } else if (format == 2) {
    /* print theta: rows = letter; cols = position in motif; 1-digit integer */
    for (i=0; i < alength; i++) {
      /* print the letter */
      fprintf(outfile, "%-*.*s%*c  ", IND, IND, 
        (i==0 ? "Simplified" : i==1 ? "pos.-specific" : i==2 ? 
	  "probability" : i==3 ? "matrix" : "" ), IND2, unhash(i));	
      for (j=0; j < w; j++) {
	int k = NINT(10.0 * theta(j,i));	/* round to 1 digit */
        if (k == 0) {
          fprintf(outfile, ":");		/* print 0 as colon */
        } else {
          fprintf(outfile, "%1x", k);		/* print 1 digit */
        }

      }
      fprintf(outfile, "\n");
    }
  }

  fprintf(outfile, "\n");
}

/**********************************************************************/
/*
	print_zij
*/
/**********************************************************************/
extern void print_zij(
  DATASET *dataset,			/* the dataset */
  MODEL *model				/* the model */
)
{
  int i, j, k;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  FILE *out=stdout;

  fprintf(out, "z_ij: lambda=%f ll=%f\n", model->lambda, model->ll);
  for (i=0; i<n_samples; i++) {			/* sequence */
    int lseq = samples[i]->length;
    int w = model->w;
    fprintf(out, ">%s\nz : ", samples[i]->sample_name);
    for (j=0; j<lseq-w+1; j++) {		/* position */
      int zij = NINT(10 * samples[i]->z[j]);	/* round z */
      fprintf(out, "%1x", zij);
      /*fprintf(out, "%11.8g ", samples[i]->z[j]);*/
    }
    fprintf(out, "\n");
    /* print sz if more than one strand */
    if (model->invcomp) {
      for (k=0; k<2; k++) {
        printf("s%d: ", k);
	for (j=0; j<lseq-w+1; j++) {		/* position */
	  int sijk = NINT(10 * samples[i]->sz[k][j]);	/* round sz */
	  printf("%1x", sijk);
	} /* site */
        printf("\n");
      } /* strand */
    }
  } /* sequence */
  printf("\n");
} /* print_zij */

/**********************************************************************/
/*
	print_wij
*/
/**********************************************************************/
extern void print_wij(
  DATASET *dataset			/* the dataset */
)
{
  int i,j;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  printf("w_ij:\n");
  for (i=0; i<n_samples; i++) {               /* sequence */
    int len = samples[i]->length;
    double *weights = samples[i]->weights;
    printf(">%s\n", samples[i]->sample_name);
    for (j=0; j<len; j++) {                   /* position */
      int w = NINT(10 * weights[j]);
      printf("%1x", w);
    }
    printf("\n");
  }
  printf("\n");
}

/**********************************************************************/
/*      get_consensus

        Get the consensus string from a motif.
	For each position, N consensus letters are found.
	If no letter has 
	probability > min_prob, 'x' is written for the first consensus
	letter and ' ' in the others.
        Otherwise, N letters are written in decreasing order of
	probability until one with min_prob is reached, and then ' ' is 
	written for the remaining letters.
*/
/**********************************************************************/
extern char *get_consensus(
  THETA theta,			/* motif theta */
  int w,			/* width of motif */
  DATASET *dataset,		/* the dataset */
  int N,			/* number of letters for each position */
  double min_prob		/* minimum cumulative prob for N letters */
)
{
  int i, j, n;
  int alength = dataset->alength;
  char *alphabet = dataset->alphabet;
  char *string = NULL;

  Resize(string, w*N+2, char);

  for (i=0; i < w; i++) {		/* position in motif */
    int maxj[MAXDEPTH];			/* array of max indices in Theta */

    /* find N letters at position i with highest probability (in order) */
    for (n = 0; n < N; n++) {		/* current depth */	
      double max = LITTLE;		/* current max probability */
      for (j=0; j < alength; j++) {	/* letter */
	if (theta(i, j) > max) {
	  max = theta(i, j);		/* maximum probability */
	  maxj[n] = j;			/* current letter with n-th best prob */
	}
      }
      theta(i, maxj[n]) = -theta(i, maxj[n]);	/* tag this position as used */
    } 

    /* restore theta */
    for (n = 0; n < N; n++) {		/* current depth */	
      theta(i, maxj[n]) = -theta(i, maxj[n]);   /* untag */	
    }

    /* set up the consensus strings for position i */
    for (n = 0; n < N; n++) { 			/* current depth */	
      if (theta(i, maxj[n]) < min_prob) {
        string[(n*w)+i] = (n==0 ? 'x' : ' ');	/* below cutoff */
      } else { 
        string[(n*w)+i] = alphabet[maxj[n]]; 	/* set n'th consensus */
      }
    }
  }
  string[((N-1)*w)+i] = '\0';			/* terminate string */

  return string;
} /* get_consensus */

/**********************************************************************/
/*
	print_candidates

	Print the candidate models found.

*/
/**********************************************************************/
static void print_candidates(
  CANDIDATE *candidates, 	/* list of candidate models */
  DATASET *dataset,		/* the dataset */
  int max_w			/* maximum width for motifs */
)
{
  int i, w;
  int hdr = 35;
  int tail = PAGEWIDTH - hdr;
  double m, e;			/* for printing signficance */

  printf("\nCandidate motifs:\n");
  printf("width nsites  ll        signif     consensus\n");
  printf("----- ------  --------- ---------- ---------\n");

  for (w=1; w<=max_w; w++) {

    if (candidates[w].sig > 1) continue;	/* skip unused candidates */

    exp10_logx(candidates[w].sig/log(10.0), m, e, 3);

    printf ("%5d %6.1f %c%9.0f %5.3fe%+04.0f ", 
      w, 
      candidates[w].lambda * wps(dataset, w), 
      (candidates[w].pal ? 'P' : ' '),
      candidates[w].ll,
      m, e);
    printf("%-*.*s\n", tail, tail, candidates[w].cons);
    for (i=tail; i < w; i+=tail) printf("%*.*s%-*.*s\n", 
      hdr, hdr, "", tail, tail, candidates[w].cons+i);
  }
} /* print_candidates */

/**********************************************************************/
/*
	print_dataset_summary
*/
/**********************************************************************/
extern void print_dataset_summary (
  DATASET *dataset 			/* structure containing sequences */
)
{
  int i, pcol;
  int w = MSN + 15;
  char *datafile = dataset->datafile;	/* name of the training set file */
  char *alphabet = dataset->alphabet;	/* alphabet of sequences */
  char *negfile	= dataset->negfile;	/* name of negative example file */

  /* set up printing spacers */
  Resize(stars, PAGEWIDTH+1, char);
  for (i=0; i<PAGEWIDTH; i++) {
    stars[i] = '*';
  }
  stars[PAGEWIDTH] = '\0';

  /* announce the training set */
  printf("%s\nTRAINING SET\n%s\n", stars, stars);

  /* print name of file and alphabet */
  printf(
"DATAFILE= %s\n"
"ALPHABET= %s\n", datafile, alphabet
  );

  /* print name of negative dataset */
  if (negfile){ printf("NEGATIVES= %s\n", negfile); }

  /*
    print a table of sequence lengths 
  */

  /*   print table header */
  for (pcol = w; pcol < 80; pcol += w) {
    printf("%-*.*s %6s %6s%2s", MSN, MSN, 
      "Sequence name", "Weight", "Length", " ");
  }
  printf("\n");
  for (pcol = w; pcol < 80; pcol += w) {
    printf("%-*.*s %6s %6s%2s", MSN, MSN, 
      "-------------", "------", "------", " ");
  }
  printf("\n");

  /*   print table columns */
  pcol = 0; 
  for (i=0; i<dataset->n_samples; i++) {
    SAMPLE *sample = dataset->samples[i];
    char *sample_name = sample->sample_name;
    double wgt = sample->sw;
    int lseq = sample->length;
    /* print the sample name and its length */
    pcol += w;          /* new line for print out? */
    if (pcol >= 80) { printf("\n"); pcol = w; }
    printf("%-*.*s %6.4f %6d%2s", MSN, MSN, sample_name, wgt, lseq, " ");
  }

  /* finish section */
  printf("\n%s\n\n", stars);
} /* print_dataset_summary */

/**********************************************************************/
/*
	print_command_summary

	Print the command line summary
*/
/**********************************************************************/
extern void print_command_summary(
  MODEL *model,			/* the model */
  DATASET *dataset 		/* the dataset */
)
{
  int i, pcol;
  char evt_string[12];

  if (dataset->evt == BIG) {
    strcpy(evt_string, "inf");
  } else {
    sprintf(evt_string, "%8g", dataset->evt);
  }

  fprintf(stdout, "%s\nCOMMAND LINE SUMMARY\n%s\n"
"This information can also be useful in the event you wish to report a\n"
"problem with the MEME software.\n\n"
"command: %s\n\n"
"model:  mod=      %8s    nmotifs=  %8d    evt=      %8s\n"
"object function=  %s\n"
"width:  minw=     %8d    maxw=     %8d    minic=    %8.2f\n",
    stars, stars,
    dataset->command, 
    dataset->mod, dataset->nmotifs, evt_string,
    (dataset->objfun == Pv ? 
       "P-value of product of p-values" : 
       "E-value of product of p-values"),
    model->min_w, model->max_w, dataset->min_ic);
  if (dataset->ma_adj) {
    fprintf(stdout,
"width:  wg=       %8g    ws=       %8g    endgaps=  %8s\n",
      dataset->wg, dataset->ws, yesno[dataset->endgaps]);
  }
  fprintf(stdout,
"nsites: minsites= %8g    maxsites= %8g    wnsites=  %8g\n"
"theta:  prob=     %8g    spmap=    %8s    spfuzz=   %8g\n" 
"em:     prior=   %9s    b=        %8g    maxiter=  %8d\n"
"        distance= %8g\n"
"data:   n=        %8d    N=        %8d\n", 
    model->min_nsites, model->max_nsites, dataset->wnsites,
    dataset->prob, dataset->mapname, dataset->map_scale,
    dataset->priorname, dataset->beta, dataset->maxiter,
    dataset->distance,
    dataset->total_res, dataset->n_samples
  );
  if (!strcmp(dataset->alphabet, DNA0)) {
    fprintf(stdout, "strands: +");
    if (model->invcomp) fprintf(stdout, " -");
  }
  fprintf(stdout,
"\n"
"sample: seed=     %8d    seqfrac=  %8g\n", dataset->seed, dataset->seqfrac);
  if (dataset->plib_name) {
    fprintf(stdout, "Dirichlet mixture priors file: %s\n", dataset->plib_name);
  }
  /* print dataset frequencies of letters in alphabet */
  fprintf(stdout, "Letter frequencies in dataset:\n");
  for (i=0, pcol=0; i<dataset->alength; i++) {
    pcol += 8;          /* start of next printed thing */
    if (pcol >= PAGEWIDTH) {pcol=8; fprintf(stdout, "\n");}
    fprintf(stdout, "%c %5.3f ", dataset->alphabet[i], dataset->res_freq[i]);
  }
  /* print background frequencies of letters in alphabet */
  fprintf(stdout, "\nBackground letter frequencies (from %s):\n",
    dataset->bfile ? dataset->bfile : "dataset with add-one prior applied");
  for (i=0, pcol=0; i<dataset->alength; i++) {
    pcol += 8;          /* start of next printed thing */
    if (pcol >= PAGEWIDTH) {pcol=8; fprintf(stdout, "\n");}
    fprintf(stdout, "%c %5.3f ", dataset->alphabet[i], dataset->back[i]);
  }
  fprintf(stdout, "\n%s\n", stars);

} /* print_command_summary */

/**********************************************************************/
/*
	print_meme_doc

	Print documentation for MEME.

	Note: file meme.doc is taken by cut-and-paste from 
	meme-output.html.
*/
/**********************************************************************/
extern void print_meme_doc(
  char *meme_directory			/* meme source directory */
)
{
  char *name = NULL;
  int s1;
  FILE *doc;

  /* make the name of the meme documentation file */
  s1 = strlen(meme_directory);
  Resize(name, s1 + 10, char);
  strcpy(name, meme_directory);
  strcpy(name + s1, "/etc/meme.doc");

  /* open the documentation file */
  doc = fopen(name, "r");
  if (doc == NULL) {
    fprintf(stderr, "Unable to open MEME's documentation file `%s'.\n", name);
    exit(1);
  }

  printf("\n");
  PSTARS;
  printf("EXPLANATION OF RESULTS\n");
  PSTARS;

  /* copy the documentation file to standard output */
  cpyfile(doc, stdout);

  PSTARS;
} /* print_meme_doc */

/**********************************************************************/
/*
	meme_score_sequence

	Compute the sequence score for a motif.

	Returns the sequence score.
*/
/**********************************************************************/
static double meme_score_sequence(
  char *eseq,		/* integer-coded sequence to score */
  int length,		/* length of the sequence */
  int w, 		/* width of motif */
  LOGODDS logodds1, 	/* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2  	/* log-odds matrix: LOG2(m_ij/n_ij) */
)
{
  int i, j;
  double best = LITTLE;			/* sequence score */
  double score, sc1, sc2;
  double loge2 = log(2);

  /* score the sequence with motif */
  for (i=0; i <= length - w; i++) {	/* site start */
    /* calculate score of subsequence */
    for (j=0, sc1=0, sc2=0; j<w; j++) {	/* position in sequence */
      sc1 += logodds1(j, (int) eseq[i+j]);
      if (logodds2) sc2 += logodds2(j, (int) eseq[i+j]);
    } /* subsequence */
    score = logodds2 ? -LOGL_SUM(-sc1*loge2, -sc2*loge2)/loge2 : sc1;
    best = MAX(score, best);
  } /* site start */

  return best;

} /* meme_score_sequence */

/**********************************************************************/
/*
	get_thresh

	Get the optimal threshold for minimizing classification error
	by classifying positive and negative data using the motif, 
	sorting and finding the minimum error.

	Returns optimal threshold, error rate and ROC.
*/
/**********************************************************************/
static ACCURACY *get_thresh(
  int w,				/* width of motif */
  LOGODDS logodds1,			/* log-odds matrix: LOG2(m_ij/b_j) */
  LOGODDS logodds2,			/* log-odds matrix: LOG2(m_ij/a_ij) */
  DATASET *pos,				/* positive examples */
  DATASET *neg,				/* negative examples */
  BOOLEAN print_scores			/* print sorted scores */
)
{
  int i, class, iseq;
  int err;				/* sum of false pos. and false neg. */
  int min_pos, max_pos;			/* best cutoff index in sorted list */
  int best_err;				/* best classification error */
  double thresh;			/* best threshold */
  SORTED_SCORE *scores=NULL;		/* array of class/score */
  int npos = pos->n_samples;
  int nneg = neg->n_samples;
  int nseqs = npos + nneg;		/* number of sequences */
  DATASET *dataset;
  /* for ROC */
  double roc;				/* receiver operating characteristic */
  double tpp, fpp;			/* true, false positive proportions */
  double tp, fp;			/* true, false positives so far */
  double newtpp, newfpp;
  ACCURACY *acc = NULL;
  double minposscore;			/* minimum score of a positive */ 
  double maxnegscore;			/* maximum score of a negative */ 

  /* allocate space for accuracy record */
  Resize(acc, 1, ACCURACY);
 
  /* allocate space for scores */
  Resize(scores, nseqs, SORTED_SCORE);

  /* score sequences */ 
  for (class=0, iseq = 0; class<2; class++) {
    if (class) dataset = pos; else dataset = neg;
    for (i=0; i<dataset->n_samples; i++) {
      SAMPLE *s = dataset->samples[i];	/* sample */
      char *eseq = s->res;			/* integer-coded sequence */
      int lseq = s->length;			/* length of sequence */
      if (lseq < w) continue;			/* sequence too short */
      scores[iseq].class = class;
      scores[iseq].score = meme_score_sequence(eseq, lseq, w, logodds1, 
        logodds2);
      scores[iseq].id = s->sample_name;
      iseq++;
    }
  }

  /* sort sequences by score in descending order */
  qsort(scores, nseqs, sizeof(SORTED_SCORE), s_compare);

  /* 
    get threshold with minimum classification error 
    If a range of thresholds give the same error, choose the average.
  */
  roc = tpp = fpp = tp = fp = 0 ;	/* for ROC */
  best_err = err = pos->n_samples;	/* all false negatives */
  min_pos = 0;				/* threshold must be below this */
  max_pos = 0;				/* threshold must be above this */
  minposscore = BIG;			/* smallest score of positives */
  maxnegscore = -BIG;			/* smallest score of negatives */
  for (i=0; i<nseqs; i++) {
    if (scores[i].class) {			/* positive */
      err--;					/* one fewer fn */
      tp++;
      minposscore = scores[i].score;
    } else {					/* negative */
      err++;					/* one more fp */
      fp++;
      maxnegscore = MAX(maxnegscore, scores[i].score);
    }
    if (err < best_err) { 			/* start new range */
      best_err = err;
      min_pos = max_pos = i;
    } else if (err == best_err) {		/* extend current range */
      max_pos = i;
    }
    /* ROC trapezoidal rule : (y2 - y1)/2 dx */
    newtpp = tp / npos;
    newfpp = fp / nneg;
    roc += .5 * (newtpp + tpp) * (newfpp - fpp);
    tpp = newtpp;
    fpp = newfpp;
  }
  max_pos = MIN(max_pos+1, nseqs-1);
  thresh = (scores[min_pos].score + scores[max_pos].score)/2;

  /* normalize by fpp to get ROC */
  if (fpp == 0) {
    roc = 1.0;
  } else {
    roc /= fpp;
  }
  
  /* add difference between positives and negatives if ROC is 1.0 */
  if (roc == 1.0) roc += minposscore - maxnegscore;

  /* print the sorted list */
  if (print_scores) {
    printf("ROC= %f\n", roc);
    for (i=0; i<nseqs; i++) printf("%-*.*s %1d %g\n", 
      MSN, MSN, scores[i].id, scores[i].class, scores[i].score); 
  }
  
  acc->thresh = thresh;
  acc->err = best_err;
  acc->roc = roc;

  return acc;
} /* get_thresh */

/**********************************************************************/
/*
        s_compare
 
        Compare two scores in descending order.  Return <0 >0
        if the first is <, > the second.  If they are equal,
        resolves ties by returning <0 if the first has smaller class.
*/
/**********************************************************************/
static int s_compare(
  const void *v1,
  const void *v2
)
{
  const SORTED_SCORE * s1 = (const SORTED_SCORE *) v1;
  const SORTED_SCORE * s2 = (const SORTED_SCORE *) v2;
  double diff = s1->score - s2->score;
  if (diff == 0) diff = (double) (s1->class - s2->class);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* s_compare */


/**********************************************************************/
/*
	get_q

	Get the value of q which gives the optimal ROC on the training
	set.
*/
/**********************************************************************/
static double get_q(
  int nsteps,					/* try nsteps from 0 to 1 */
  int window,					/* smoothing window radius */
  int w,					/* width of motif */
  THETA theta,					/* motif theta */
  THETA neg_theta, 				/* anti-motif theta */
  double *back,					/* background motif */
  DATASET *dataset,				/* the dataset */
  DATASET *neg_dataset,				/* negative examples */
  char *str_space 		/* space for printing strand direction */
)
{
  int i, j;
  double *roc = NULL;				/* array to hold roc(q) */
  double smooth_roc;				/* smoothed value of roc */
  double best_roc;				/* maximum roc */
  double q=0;					/* mixing parameter */
  double incr = 1.0/nsteps;			/* increment for q */
  LOGODDS logodds;				/* motif log-odds matrix */
  ACCURACY *acc;				/* get_thresh struct */
  int alength = dataset->alength;		/* length of alphabet */

  /* create ROC array */
  Resize(roc, nsteps+1, double);

  /* 
    get roc for different values of q 
  */
  for (i=0; i<=nsteps; i++) {
    q = i * incr;
    logodds = make_log_odds(theta, neg_theta, back, q, w, alength);
    acc = get_thresh(w, logodds, NULL, dataset, neg_dataset, FALSE);
    roc[i] = acc->roc;				/* save roc for this q */
    myfree(acc);				/* free up space */
  } /* get roc */

  /* 
    smooth ROC and find q that gives maximum
  */
  best_roc = 0;
  for (i=0; i<=nsteps; i++) {
    double avg = 0;
    int cnt = 0;
    for (j=MAX(0,i-window); j<=MIN(nsteps, i+window); j++) {
      avg += roc[j];
      cnt++;
    }
    smooth_roc = avg/cnt;
    if (smooth_roc > best_roc) {
      best_roc = smooth_roc;
      q = i * incr;
    }
    /*printf("q= %8.3f smoothed_roc= %8.5f\n", i*incr, smooth_roc);*/
  } /* smooth ROC and get max */

  myfree(roc);					/* release space */

  printf("Q= %8.3f ROC= %8.3f\n", q, best_roc); 

  return q;
} 

/**********************************************************************/
/*
	print_sites

	Print the sites making up the model.
*/
/**********************************************************************/
static void print_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  int format,				/* 0=BLOCKS; 1=FASTA */
  char *com,				/* comment to append */
  FILE *outfile				/* where to print */
)
{
  int i, j;
  int w = model->w;			/* width of motif */
  P_PROB sites = model->maxima;		/* sites "defining" model */
  int n = model->nsites_dis;		/* number of sites */
  char *ftype = (format==0 ? "BLOCKS" : "FASTA");

  /* print header */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d in %s format%s\n", model->imotif, ftype, com);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);

  /* print the sites */
  if (format == 0) fprintf(outfile, "BL   MOTIF %d width=%d seqs=%d\n", 
    model->imotif, w, n);
  for (i=0; i<n; i++) {			/* print each site */
    int seqno = sites[i].x;		/* sequence number */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    BOOLEAN ic = sites[i].ic;		/* strand direction */
    int y = sites[i].y;			/* location of site */
    int off = ic ? s->length-w-y : y;	/* - strand offset from rgt. */
    char *res = ic ? s->resic+off : s->res+off;       /* integer sequence */
    double weight = s->sw;		/* sequence weight */
    /*double weight = sites[i].prob;*/

    /* print sequence name and position of site */
    if (format == 0) {			/* BLOCKS format */
      fprintf(outfile, "%-*.*s ( %4d) ", MSN, MSN, s->sample_name, y+1);
    } else {				/* FASTA format */
      fprintf(outfile,">%-*.*s pos %4d\n", MSN, MSN, s->sample_name, y+1);
    }

    /* print site */
    for (j=0; j<w; j++) { fputc(unhash(res[j]), outfile); }
    if (format == 0) {			/* BLOCKS format */
      fprintf(outfile, "  %g ", weight);
    }
    fputc('\n', outfile);
  } /* print each site */
  if (format == 0) {
    fprintf(outfile, "//\n\n");
  }
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile);
  fprintf(outfile, "\n\n");

} /* print_sites */

/**********************************************************************/
/*
	print_summary

	Print the summary of all the motifs.
*/
/**********************************************************************/
extern void print_summary(
  MODEL *model,				/* the model */
  DATASET *dataset,			/* the dataset */
  LO **los,				/* the LO structures */
  int nmotifs,				/* number of motifs */
  double **pv,				/* p-value of score distribution */
  FILE *outfile				/* where to print */
)
{
  /* print the motif block diagrams using all the motifs */
  fprintf(outfile, "\n\n%s\nSUMMARY OF MOTIFS\n%s\n\n", stars, stars);
  print_block_diagrams(model, dataset, los, nmotifs, pv, outfile);
  fprintf(outfile, "%s\n\n", stars);
} /* print_summary */

/**********************************************************************/
/*
	score_sites

	Score and get the pvalues of the sites in the model.
	Sort in order of increasing p-value.
*/
/**********************************************************************/
static void score_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  LO *lo,				/* LO structure */
  double *pv 				/* p-values for scores of this motif */
)
{
  int isite;
  P_PROB sites = model->maxima;		/* sites "defining" model */
  int n = model->nsites_dis;		/* number of sites */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  STYPE stype = invcomp ? Combine : Protein;	/* Protein works for DNA too */
  SCORE **scores = NULL;		/* the site scores */
  int old_seqno = -1;

  for (isite=0; isite<n; isite++) {	/* site */
    int seqno = sites[isite].x;		/* sequence number */
    int y = sites[isite].y;		/* location of site */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;		/* length of sequence */
    char *seq = s->seq;			/* the ascii sequence */
    double pvalue; 			/* score p-value */
 
    /* score the sequence if new */
    if (old_seqno != seqno) {
      BOOLEAN xlate_dna = FALSE;	/* not xlating */
      if (old_seqno >= 0) free_2array(scores, 1);
      scores = score_sequence(stype, xlate_dna, seq, lseq, 1, &lo);
      old_seqno = seqno;
      s->minpv = 1.0;
    }

    pvalue = pv[(int) scores[0][y].score];	/* p-value */

    /* save MINUS the p-value in the .prob field of sites */
    sites[isite].prob = -pvalue;

    /* update minimum p-value of sites */
    if (pvalue < s->minpv) s->minpv = pvalue;

  } /* get p-values of sites */
  free_2array(scores, 1);               /* free space */

  /*
    sort the sites by p-value
  */
  qsort((char *) sites, n, sizeof(p_prob), pY_compare);

  /*
    change sign of p-values back
  */
  for (isite=0; isite<n; isite++) sites[isite].prob *= -1;

} /* score_sites */

/**********************************************************************/
/*
	print_site_diagrams

	Make block diagrams of the actual sites in the model
	and print them.
	Sequences are sorted by the minimum p-value of sites in them.
*/
/**********************************************************************/
static void print_site_diagrams(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  int nmotifs,				/* number of motifs in los */
  LO *los[MAXG],			/* logodds structure for motif */
  FILE *outfile				/* where to print */
)
{
  int i, j, isite;
  P_PROB sites = model->maxima;		/* sites "defining" model */
  int n = model->nsites_dis;		/* number of sites */
  int nseqs = dataset->n_samples;	/* number of sequences in dataset */
  BOOLEAN dna = dataset->dna;		/* dataset is DNA if true */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  BOOLEAN xlate_dna = FALSE;		/* don't translate */
  BOOLEAN best_motifs = FALSE;		/* use all sites */
  double m_thresh = 1;			/* show all sites as strong */
  STYPE stype = dna ? (invcomp ? Combine : Norc) : Protein;
  int nseqs_with_sites;			/* number of sequences with sites */
  int *seqlist = NULL;			/* ordered list of sequences w/sites */
  int *hits = NULL;			/* store hits */
  double *pvalues = NULL;		/* store pvalues */
  char *f = "%-*s%s %8s  %s\n";		/* format */

  /*
    Create the list to contain sequences sorted by minimum p-value
  */
  Resize(seqlist, nseqs, int);

  /*
    Clear list of sites for each sequence 
  */
  for (i=0; i<nseqs; i++) dataset->samples[i]->nsites = 0;

  /*
    Find which sequences have sites and create list of sites for each
  */
  for (isite=nseqs_with_sites=0; isite<n; isite++) {	/* site */
    int seqno = sites[isite].x;		/* sequence number */
    int y = sites[isite].y + 1;		/* location of site (plus 1) */
    int ic = sites[isite].ic;		/* site on reverse complement strand */
    SAMPLE *s = dataset->samples[seqno];/* sequence */

    /* record the sequence as containing a site */
    if (!s->nsites) seqlist[nseqs_with_sites++] = seqno;

    /* store the site in its list of sites */
    Resize(s->sites, s->nsites+1, int);
    s->sites[(s->nsites)++] = ic ? -y : y;	/* +/- site offset by 1 */
  } /* site */

  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "\tMotif %d block diagrams\n", nmotifs);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, f, MSN, "SEQUENCE NAME", "", "POSITION P-VALUE", "MOTIF DIAGRAM");
  fprintf(outfile, f, MSN, "-------------", "", "----------------", "-------------");

  /*
    create and print a block diagram for each sequence with sites 
  */
  for (i=0; i<nseqs_with_sites; i++) {	/* sequence */
    int seqno = seqlist[i];		/* current sequence */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;		/* length of sequence */
    char *name = s->sample_name;	/* name of sequence */
    double minpv = s->minpv;		/* minimum p-value of sites */
    char hdr[MSN+20];			/* header */
    TILING tiling;			/* tiling struct */

    /*
      create storage for hits and pvalues and clear them
    */
    Resize(hits, lseq, int);
    Resize(pvalues, lseq, double);
    for (j=0; j<lseq; j++) { hits[j] = 0; pvalues[j] = 0; }

    /* copy hits from s->nsites into hits array */
    for (j=0; j<s->nsites; j++) {
      int y = abs(s->sites[j]) - 1;	/* position of site */
      int m = (s->sites[j] > 0) ? los[nmotifs-1]->name : -los[nmotifs-1]->name;
      hits[y] = m;			/* +/- motif */
    }

    /* put the hits in TILING struct */
    tiling.hits = hits;
    tiling.pvalues = pvalues;

    /* create the block diagram */
    tiling.diagram = create_diagram(dna, stype, xlate_dna, best_motifs, FALSE,
      m_thresh, nmotifs, los, lseq, FALSE, tiling);

    /* print the diagram */
    sprintf(hdr, "%-*.*s %16.2g  ", MSN, MSN, name, minpv);
    print_diagram(tiling.diagram, hdr, outfile);

    myfree(tiling.diagram);		/* release space */
  } /* sequence */

  /* print a final line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");
  
  myfree(seqlist);
  myfree(hits);
  myfree(pvalues);

} /* print_site_diagrams */

/**********************************************************************/
/*
	align_sites

	Align all sites that make up the model.
*/
/**********************************************************************/
static void align_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  LO *lo,				/* LO structure */
  double *pv,				/* pvalues for scores of this motif */
  FILE *outfile				/* stream for output */
)
{
  int i, ii, isite;
  int w = model->w;			/* length of site */
  P_PROB sites = model->maxima;		/* sites "defining" model */
  int n = model->nsites_dis;		/* number of sites */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  int imotif = lo->name;		/* name of motif */
  char site[MAXSITE+1], pre[10+1], post[10+1];
  /* added by M.H. */	
  double aveSS = 0;
  int numHits = 0;

  /* print header */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile,
    "\tMotif %d sites sorted by position p-value\n", imotif);
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fputc('\n', outfile);
  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "Sequence name", 
    invcomp ? " Strand" : "");

  /* added by M.H. */
  /* remove \n and add a > 0 ? statement */
  fprintf(outfile, "%6s %9s %10s %*sSite%*s", 
    "Start", "P-value", "", ( (w/2-2) >0 ? (w/2-2) : 0), "", ((w - w/2 - 4) > 0 ? (w - w/2 - 4) : 0), "");
/*  fprintf(outfile, "%6s %9s %10s %*sSite%*s\n", 
    "Start", "P-value", "", w/2 - 2, "", w - w/2 - 4, "");*/
  /* added by M.H. */
  /* output 'SecStruc' */
  if (dataset->secondaryStructureFilename != NULL) {
	 	fprintf(outfile, "%*s %10s %8s\n", MAX(0, w - ( (w-w/2-4) > 0 ? (w-w/2-4) : 0 ) - ( (w/2-2) >0 ? (w/2-2) : 0) - 4 ) ,  "","","SecStruc");
	 }else{
	 	fprintf(outfile, "\n");
	}

  fprintf(outfile, "%-*.*s%s ", MSN, MSN, "-------------", 
    invcomp ? " ------" : "");
  fprintf(outfile, "%6s %8s %10s ", "-----", "---------", "");
  for (i=0; i<w; i++) fputc('-', outfile);
  /* added by M.H. */
  if (dataset->secondaryStructureFilename != NULL) {
	  fprintf(outfile, " %10s %8s", "", "--------");
  }
  fputc('\n', outfile);

  /* 
    print sites that make up the model
  */
  for (isite=0; isite<n; isite++) {	/* site */
    int seqno = sites[isite].x;		/* sequence number */
    int y = sites[isite].y;		/* location of site */
    BOOLEAN ic = sites[isite].ic;	/* strand direction */
    double pvalue = sites[isite].prob;	/* position p-value */
    SAMPLE *s = dataset->samples[seqno];/* sequence */
    int lseq = s->length;		/* length of sequence */
    char *seq = s->seq;			/* the ascii sequence */
    char *sample_name = s->sample_name;	/* name of sample */

    /* print name and strand */
    fprintf(outfile, "%-*.*s%s ", MSN, MSN, sample_name, 
      invcomp ? (ic ? "     -" : "     +") : "");

    /* print position and position p-value */
    fprintf(outfile, "%6d %9.2e", y+1, pvalue);

    /* get the aligned sequence parts */
    if (!ic) {				/* + strand */
      /* pre */
      for (i=y-10, ii=0; i<y; i++) {
	if (i<0) continue;
	pre[ii++] = seq[i];
      }
      pre[ii] = '\0';
      /* site */
      for (i=y, ii=0; ii<w; i++)  site[ii++] = seq[i]; 
      site[ii] = '\0';
      /* post */
      for (i=y+w, ii=0; ii<10 && i<lseq; i++) post[ii++] = seq[i];
      post[ii] = 0;

    } else {				/* - strand */
      /* pre */
      for (i=y+w+9, ii=0; i>=y+w; i--) {
	if (i>=lseq) continue;
	pre[ii++] = comp_dna(seq[i]);
      }
      pre[ii] = '\0';
      /* site */
      for (i=y+w-1, ii=0; ii<w; i--) site[ii++] = comp_dna(seq[i]);
      site[ii] = '\0';
      /* post */
      for (i=y-1, ii=0; ii<10 && i>=0; i--) post[ii++] = comp_dna(seq[i]);
      post[ii] = '\0';
    } /* strand */

	 /* added by M.H. */
	 /* remove \n  */
    /* print the alignment */
    if (pre[0] == '\0') {			/* print a dot in empty pre */
      fprintf(outfile, " %10s %-*s %-10s", ".", w, site, post);
    } else {
      fprintf(outfile, " %10s %-*s %-10s", pre, w, site, post);
    }
/*    if (pre[0] == '\0') {			/* print a dot in empty pre *
      fprintf(outfile, " %10s %-*s %-10s\n", ".", w, site, post);
    } else {
      fprintf(outfile, " %10s %-*s %-10s\n", pre, w, site, post);
    } */

	/* added by M.H.. */
	/* print EF/PU value and compute average */
    if (dataset->secondaryStructureFilename != NULL) {
	 	fprintf(outfile, " %f\n",s->ss_value[y]);
	 	aveSS += s->ss_value[y];
		numHits ++;
	 }else{
	 	fprintf(outfile, "\n");
	 }


  } /* site */

  /* added by M.H.. */
  /* print average */
  if (dataset->secondaryStructureFilename != NULL) {
  		aveSS /= numHits;
	 	fprintf(outfile, "average single-strandedness: %f\n",aveSS);
  }

  /* print line of hyphens */
  for (i=0; i<PAGEWIDTH; i++) fputc('-', outfile); fprintf(outfile, "\n\n");

} /* align_sites */

