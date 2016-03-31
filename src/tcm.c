/*
 * $Id: tcm.c,v 1.2 2005/10/25 19:06:39 nadya Exp $
 * 
 * $Log: tcm.c,v $
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:27:46  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 9-30-99 tlb; remove computation of cum. background probability */
/* 7-02-99 tlb; remove clobbering of theta */
/* 6-23-99 tlb; rewrite to support alternate DNA strands */
/*	
	EM algorithm.

	Two component mixture model. 
*/
	
#include "meme.h"

static double smooth(
  int w,				/* width to smooth over */
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
);

/**********************************************************************/
/*
	tcm_e_step	

	Do the E step of EM.

	Estimate the expectation of model 1 for each position in the data.
	In other words, calculate E[z_ij] in z.

	Updates z.  

	Returns log pr(X | theta, lambda).

	Time: O(n_samples*lseq*w)
*/
/**********************************************************************/
double tcm_e_step(
  MODEL *model,			/* the model */
  DATASET *dataset  		/* the dataset */
)
{
  int i, j, k, ii;
  THETA logtheta1 = model->logtheta;	/* motif log(theta) */
  int w = model->w;			/* motif width */
  int n_samples = dataset->n_samples;	/* number of sequences */
  BOOLEAN invcomp = model->invcomp;     /* use reverse complement strand, too */
  int ndir = invcomp ? 2 : 1;           /* number of strands */
  double log_sigma = log(1.0/ndir);	/* log \sigma */
  double lambda = model->lambda;	/* \lambda of tcm model */
  double log_lambda = LOG(lambda);	/* log \lambda */
  double log_1mlambda = LOG(1-lambda);	/* log (1 - \lambda) */
  double logpX;				/* log likelihood; no erase or smooth */

  /* E step */

  convert_theta_to_log(model, dataset);

  /* calculate all the posterior offset probabilities */
  logpX = 0;

  for (i=0; i < n_samples; i++) {	/* sequence */
    SAMPLE *s = dataset->samples[i];
    int lseq = s->length;
    int m = lseq - w + 1;		/* number of possible sites */
    double *zi = s->z;			/* Pr(z_ij=1 | X_i, \theta) */
    double *not_o = s->not_o;		/* Pr(v_ij = 1) */
    double *lcb = s->logcumback;	/* cumulative background probability */

    if (lseq < w) continue;		/* sequence too short for motif */

	 /* added by M.H. */
	 /* use log sigma_ij * lambda * m instead of lambda if secondary structure information is given */
	 /* NOTE: log_sigma is the prior for + or - strand --> here only + strand --> log_sigma = 0 */
	 if (dataset->secondaryStructureFilename != NULL) {
	 
	 	/* first check if the maximum sigma * (lambda*m) > 1   --> if so P(Zij=1 | \phi) can be > 1 */
	   double Pcount = dataset->secondaryStructurePseudocount;
		double maxPrior = ((s->max_ss_value + Pcount) / (s->sum_ss_value + (m * Pcount))) * (lambda*m);

		if (maxPrior > 1.0 ) {
			/* compute new pseudocount that gives sigma_i_max = 1 = (max_ss_value + pseudocount) / (\sum (ss_value[i] + pseudocount)) * lambda * m */
			Pcount = (-1 * s->max_ss_value * lambda * m + s->sum_ss_value) / (m * (lambda - 1));
	
			/* for statistics keep maximum adjustment */
			if (Pcount - dataset->secondaryStructurePseudocount > MAXADJUST) {
				MAXADJUST = Pcount - dataset->secondaryStructurePseudocount;
			}
		}
	 	
		/* compute new sigmas with this pseudocount */
		double sum = s->sum_ss_value + (m * Pcount);  		 /* \sum ss_value[i] + m*pseudocount */
    	for (j=0; j < m; j++) {
	   	s->sigma[j] = (s->ss_value[j] + Pcount) / sum;
		}
	 }


    for (k=0; k<ndir; k++) {		/* strand */
      BOOLEAN ic = (k==1);		/* doing - strand */
      double *szik = s->sz[k];		/* Pr(X_i | z_ij=1, s_ijk=1, \theta) */

      for (j=0; j<m; j++) {		/* site start */
		
	 	   /* added by M.H. */
		   /* use the prior instead of lambda */
	 	   if (dataset->secondaryStructureFilename != NULL) {
			  double p = MIN(1.0, ( s->sigma[j] * lambda * m ) ) ;		/* rounding */
			  log_lambda = LOG(p);	
  			  log_1mlambda = LOG(1-p);	
		   }
		
       	 /* log Pr(X_ij | s_ijk=1, \theta0) \sigma (1-\lambda) */
			 double log_pXijtheta0 = log_sigma + log_1mlambda;	
         /* log Pr(X_ij | s_ijk=1, \theta1) \sigma \lambda */
			double log_pXijtheta1 = log_sigma + log_lambda;
         int off = ic ? lseq-w-j : j;	/* - strand offset from rgt. */
         char *res = ic ? s->resic+off : s->res+off;	/* integer sequence */

        /* calculate the probability of positions in the site under the
  	  background and foreground models
        */
        log_pXijtheta0 += Log_back(lcb, j, w);
        for (ii=0; ii<w; ii++) log_pXijtheta1 += logtheta1(ii, (int)res[ii]);
 
        /* set log szik to:
          Pr(X_i | z_ij=1, s_ijk=1, \theta) \sigma \lambda
        */
        szik[j] = log_pXijtheta1;
 
        /* set z_ij to log Pr(X_ij | \phi): (6-21-99 tlb)
          log(
	    \sigma * sum_{k=0}^{ndir-1} ( 
	      Pr(X_i|z_ij=1, s_ijk=1, \theta) \lambda +
	      Pr(X_i|z_ij=0, s_ijk=1, \theta) (1-\lambda) 
	    )
          )
        */
        zi[j] = k==0 ? LOGL_SUM(log_pXijtheta0, log_pXijtheta1) : 
          LOGL_SUM(zi[j], LOGL_SUM(log_pXijtheta0, log_pXijtheta1));
      } /* site start */
    } /* strand */

    /* compute log Pr(X | \phi) = sum_i,j log(Pr(X_ij)) */
    for (j=0; j<m; j++) {			/* site start */
      logpX += zi[j];				/* z_ij = log Pr(X_ij | \phi) */
    }

    /* sz_ijk : normalize, delog and account for erasing
      Pr(z_ij=1, s_ijk=1 | X_i, \phi) \approx
           P(z_ij=1, s_ijk=1 | X_i, \phi) P(v_ij = 1)
    */
    for (k=0; k<ndir; k++) {		/* strand */
      double *szik = s->sz[k];		/* Pr(X_i | z_ij=1, s_ijk=1, \phi) */
      for (j=0; j<m; j++) {		/* site start */
        /* note zi[j] holds Pr(X_ij|\phi) */
        szik[j] = MIN(1.0, exp(szik[j] - zi[j]) * not_o[j]);	/* roundoff */
      } /* site */
    } /* strand */

    /* z_ij : sum of sz_ijk */
    for (j=0; j<m; j++) {		/* site start */
      for (k=zi[j]=0; k<ndir; k++) {	/* strand */
        zi[j] += s->sz[k][j];
      } /* strand */
      zi[j] = MIN(1.0, zi[j]);		/* avoid roundoff errors */
    } /* site */
    for (j=m; j<lseq; j++) {		/* tail of sequence */
      zi[j] = 0;
    }

  } /* sequence */

  /* smooth so no window of size w has z_i which sum to greater than 1.0 */
  (void) smooth(w, model, dataset);

  return (logpX/log(2.0));
} /* tcm_e_step */

/***********************************************************************/
/*
  smooth

  Normalize so that no local region w wide has z_ij sum of > 1.0.
  Winner-take-all: the largest value of z_ij is never reduced.

  Returns the total expected number of sites of motif.
*/ 
/***********************************************************************/
static double smooth(
  int w,				/* width to smooth over */
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
)
{
  int i, j, k, p;
  int n_samples = dataset->n_samples;
  SAMPLE **samples = dataset->samples;
  double p_sum = 0.0;
  BOOLEAN invcomp = model->invcomp;     /* use reverse complement strand, too */
  int ndir = invcomp ? 2 : 1;           /* number of strands */

  for (i=0; i<n_samples; i++) {		/* sequence */
    int ioff;
    SAMPLE *sample = samples[i];	/* sample i */
    int lseq = sample->length;
    double *zi= sample->z;		/* z */
    int max_o = lseq - w + 1;		/* largest possible offset */

    if (lseq < w) continue;		/* sequence too short for motif */

    /* normalize adjacent windows of length w, then shift and repeat */
    for (ioff = 0; ioff < MIN(w, max_o); ioff+=2) {	/* window start */
      for (j=ioff; j<max_o; j += w) {		/* adjacent windows */
	double local_z = 0.0;
        double max_z = 0;			/* find largest z_ij */
        int max_p = 0;
        int last_p = MIN(j+w, max_o);
	for (p=j; p<last_p; p++) {		/* position */
          double z = zi[p];
	  local_z += z;				/* compute local motif z sum */
	  if (z > max_z) {		
	    max_z = z;				/* largest z in window */
	    max_p = p;				/* position with largest z */
	  }
	}
	/* normalize if necessary; leave largest z in window unchanged */
	if (local_z > 1.0) {			/* normalize */
	  double scale = (1 - max_z) / (local_z - max_z);
	  for (k=0; k<ndir; k++) {		/* strand */
	    double *szik = sample->sz[k];
	    for (p=j; p<last_p; p++) {		/* position */
	      if (p != max_p) {
		if (k==0) zi[p] *= scale;	/* normalize z */
		szik[p] *= scale;		/* normalize sz */
	      }
	    } /* position */
	  } /* strand */
	} /* normalize */
      } /* adjacent windows */
    } /* window start */

    /* calculate p_sum */
    for (j=0; j < max_o; j++) {
      p_sum += zi[j];
    }

  } /* n_samples loop */

  return p_sum;
} /* smooth */

