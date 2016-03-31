/*
 * $Id: meme_util.c,v 1.2 2005/10/25 19:06:39 nadya Exp $
 * 
 * $Log: meme_util.c,v $
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:22:39  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1995, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 5-26-00 tlb; add cons0 to copy_model */
/* 5-23-00 tlb; fix finding of sites for OOPS to reflect scaling of zij */
/* 7-29-99 tlb; add nsites and nsites_obs to copy_model */
/* 7-27-99 tlb; add rentropy to MODEL; rem. p_point, nmotifs from create_model*/
/* 7-14-99 tlb; move get_sites and make_log_odds here from display.c */
/* 7-02-99 tlb; add logtheta to model, lambda_obs to copy_model */
/* 6-28-99 tlb; remove sigma from model, add psites */

#include <meme.h>

/**********************************************************************/
/*
  copy_theta

  Copy the theta matrix to another array.

*/
/**********************************************************************/
extern void copy_theta(
  THETA s,	 		/* source */
  THETA d,			/* destination */
  int w,			/* width of motif */
  int alength			/* length of alphabet */
)
{ 
  int i, j;
  for (i = 0; i < w; i++) {		/* col in motif */
    for (j = 0; j < alength; j++) {	/* row in motif */
        theta_ref(d, i, j) = theta_ref(s, i, j);
    }
  }
}

/**********************************************************************/
/*
	get_not_o

	Calculate the probability that each possible site
     	start does not overlap a previously found site;
     	this is taken to be the minimum of all the probabilities
     	of a site starting within the possible site
*/
/**********************************************************************/
extern void get_not_o(
  DATASET *dataset,			/* the dataset */
  int w,				/* width of motif */
  BOOLEAN get_log			/* compute log_not_o if true */
)
{
  int i,j,k;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;

  for (i=0; i<n_samples; i++){		/* sequence */
    SAMPLE *s = samples[i];
    int lseq = s->length;
    double *weights = s->weights;       /* prb not in a previous site */
    double *not_o = s->not_o;		/* prb not overlapping a site */
    int *log_not_o = s->log_not_o;	/* prb not overlapping a site */

    if (lseq < w) continue;		/* sequence to short for motif */

    for (j=0; j <= lseq - w; j++) { 	/* site start */
      not_o[j] = 1.0;			/* assume not overlapping */
      for (k=j; k < j+w; k++) { 	/* position in sequence */
        if (weights[k] < not_o[j]) not_o[j] = weights[k];
      }
      if (get_log) log_not_o[j] = INT_LOG(not_o[j]);
    }
    for (j=lseq-w+1; j < lseq ; j++) { 	/* beyond possible site starts */
      not_o[j] = 1;
      if (get_log) log_not_o[j] = 0;
    }
  }
} /* get_not_o */

/**********************************************************************/
/*
	create_model

	Create a model structure.
*/
/**********************************************************************/
extern MODEL *create_model(
  MTYPE mtype,				/* type of model */
  BOOLEAN invcomp,			/* use inv comp strand  */
  int max_w,				/* maximum width of motif */
  int alength				/* length of alphabet */
)
{
  MODEL *model = (MODEL *) mymalloc(sizeof(MODEL));

  model->mtype = mtype;
  model->invcomp = invcomp;

  /* add one extra column in case model is a palindrome */
  create_2array(model->theta, double, max_w+1, alength+1);	/* X */
  create_2array(model->logtheta, double, max_w+1, alength+1);	/* X */
  create_2array(model->obs, double, max_w+1, alength);
  model->maxima = NULL;
  model->w = 0;

  model->iter = 0;
  return model;
}

/**********************************************************************/
/*
  copy_model

  Copy a model structure.
*/
/**********************************************************************/
extern void copy_model(
  MODEL *m1, 				/* source */
  MODEL *m2,				/* destination */
  int alength				/* length of alphabet */
)
{
  int i;

  m2->mtype = m1->mtype;
  m2->min_w = m1->min_w;
  m2->max_w = m1->max_w;
  m2->pw = m1->pw;
  m2->min_nsites = m1->min_nsites;
  m2->max_nsites = m1->max_nsites;
  m2->psites = m1->psites;
  if (m1->maxima) {				/* copy maxima if they exist */
    Resize(m2->maxima, m1->nsites_dis+1, p_prob);
    bcopy((char *) m1->maxima, (char *) m2->maxima, 
      m1->nsites_dis*sizeof(p_prob));
  }
  m2->pal = m1->pal;
  m2->invcomp = m1->invcomp;
  m2->imotif = m1->imotif; 
  m2->w = m1->w;
  copy_theta(m1->theta, m2->theta, m1->w, alength);
  copy_theta(m1->obs, m2->obs, m1->w, alength);
  m2->lambda = m1->lambda;
  m2->lambda_obs = m1->lambda_obs;
  m2->nsites = m1->nsites;
  m2->nsites_obs = m1->nsites_obs;
  m2->nsites_dis = m1->nsites_dis;
  strcpy(m2->cons, m1->cons); 
  strcpy(m2->cons0, m1->cons0); 
  for (i=0; i<m1->w; i++) {
    m2->rentropy[i] = m1->rentropy[i];
  }
  m2->rel = m1->rel;
  m2->ll = m1->ll;
  m2->mll_0 = m1->mll_0;
  m2->mll_1 = m1->mll_1;
  m2->logpv = m1->logpv;
  m2->logev = m1->logev; 
  m2->ID = m1->ID; 
  m2->iter = m1->iter; 
  m2->iseq = m1->iseq;
  m2->ioff = m1->ioff;
} /* copy_model */

/**********************************************************************/
/*
	get_sites

	Get the principal sites contributing to a motif.

	OOPS		position with highest z_ij > 0.0 in each sequence
	ZOOPS		position with highest z_ij > 0.5 in each sequence
	TCM		all positions with z_ij > 0.5

	Assumes that z_ij is set for current motif.
	Returns a list consisting of sites and sets n, the number of sites.
*/
/**********************************************************************/
extern SITE *get_sites(
  DATASET *dataset,			/* the dataset */
  MODEL *model,				/* the model */
  int *n,				/* number of sites found */
  int *best_site			/* index of best site in array */
)
{
  int i, j;
  int nsites = 0;			/* number of sites */
  SITE *sites = NULL;			/* list of sites */
  int n_samples = dataset->n_samples;	/* number of sequences */
  SAMPLE **samples = dataset->samples;	/* the sequences */
  MTYPE mtype = model->mtype;		/* type of model */
  int w = model->w;			/* width of motif */
  BOOLEAN invcomp = model->invcomp;	/* using invcomp strand */
  double best_z = -1;			/* best overall zij */

  if (mtype==Oops || mtype==Zoops) {	/* (Z)OOPS model */

    for (i=0; i<n_samples; i++) {	/* sequence */
      SAMPLE *sample = samples[i];	/* sample */
      int lseq = sample->length;	/* length of sequence */
      double max_zij = -1;		/* flag no z_ij found */
      int max_j = -1;			/* position of site */
      double **sz = sample->sz;		/* sz_i */

      if (lseq < w) continue;		/* sequence to short for motif */

      /* find maximum z_ij */
      for (j=0; j<lseq-w+1; j++) {	/* position */
        double zij = sample->z[j];	/* z_ij */
        if (zij > max_zij) {		/* bigger found */
          max_zij = zij;
          max_j = j;
        }
      } /* position */

      /* record site */
      if ((max_zij && mtype == Oops) || max_zij > 0.5) {	/* site found */
        if (nsites % 100 == 0) Resize(sites, nsites+101, SITE);
        sites[nsites].seqno = i;
        sites[nsites].pos = max_j;
        sites[nsites].zij = max_zij;
        /* on reverse complement strand? */
        sites[nsites].invcomp = (invcomp && sz[1][max_j]>sz[0][max_j]) ? 1 : 0;
        if (max_zij > best_z) {
          *best_site = nsites;
          best_z = max_zij;
        }
        nsites++;
      } /* site found */

    } /* sequence */

  } else {				/* TCM model */

    for (i=0; i<n_samples; i++) {	/* sequence */
      SAMPLE *sample = samples[i];	/* sample */
      int lseq = sample->length;	/* length of sequence */
      double **sz = sample->sz;		/* sz_i */

      if (lseq < w) continue;		/* sequence to short for motif */

      /* find all z_ij > 0.5 */
      for (j=0; j<lseq-w+1; j++) {	/* position */
        double zij = sample->z[j];	/* z_ij */
        if (zij > 0.5) { 		/* record site */
	  if (nsites % 100 == 0) Resize(sites, nsites+101, SITE);
	  sites[nsites].seqno = i;
	  sites[nsites].pos = j;
	  sites[nsites].zij = zij;
	  /* on reverse complement strand? */
	  sites[nsites].invcomp = (invcomp && sz[1][j]>sz[0][j]) ? 1 : 0;
	  if (zij > best_z) {
	    *best_site = nsites;
	    best_z = zij;
	  }
	  nsites++;
        } /* record site */
      } /* position */
    } /* sequence */

  } /* TCM model */

  /* return results */
  *n = nsites;
  return sites;
} /* get_sites */

/**********************************************************************/
/*	
	make_log_odds

	Compute the log-odds matrix from the motif, negative motif and
	background letter frequencies.
		l = log(p/(q*n+(1-q)b))
	where 	p = positive model,
		n = negative model,
		b = background model.

	Include 'X' in the logodds matrix.
	X is given a score equal to the weighted average score for the
	column.
*/
/**********************************************************************/
extern LOGODDS make_log_odds(
  THETA theta1,			/* positive motif theta */
  THETA theta0,			/* negative motif theta; use 0 if NULL */
  double *back,			/* background frequencies; use 0 if NULL */
  double q,			/* 0<=q<=1, mixing parameter */
  int w,			/* width of motif */
  int alength 			/* length of alphabet */
)
{
  int i, j;
  LOGODDS logodds = NULL;		/* the logodds matrix */

  /* 
    compute the log-odds matrix 
  */
  Resize(logodds, w, double *);		/* create rows */ 
  for (i=0; i<w; i++) {			/* site position */
    logodds[i] = NULL;			/* create columns */
    Resize(logodds[i], alength+1, double);	/* include 'X' */

    /* calculate log-odds for this column */
    logodds(i,alength) = 0;		/* 'X' */
    for (j=0; j<alength; j++) {		/* letter */
      double p = theta1(i,j);		/* positive motif */
      double n; 			/* negative motif */
      if (!theta0) {			/* no negative motif given */
        n = back[j];			/* background */
      } else if (!back) {		/* no background given */
        n = theta0(i,j);		/* negative motif */
      } else {				/* negative and background given */
        n = q*theta0(i,j) + (1-q)*back[j];/* blend negative & background */
      }
      if (n==0) {
        logodds(i,j) = 0;		/* letter with zero prob */
      } else {
        logodds(i,j) = LOG2(p/n);
      }
      logodds(i,alength) += n * logodds(i,j);		/* 'X' */
    } /* letter */
  } /* column */

  return logodds;
} /* make_log_odds */

