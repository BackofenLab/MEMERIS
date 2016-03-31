/*
 * $Id: oops.c,v 1.3.4.1 2006/01/24 20:44:08 nadya Exp $
 * 
 * $Log: oops.c,v $
 * Revision 1.3.4.1  2006/01/24 20:44:08  nadya
 * update copyright
 *
 * Revision 1.3  2006/01/09 02:47:09  tbailey
 * Fixed bug in e_step: put parentheses around conditional in
 * 	double log_pXijtheta = LOG(not_o[j]) + (conditional)
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:23:52  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994-2006, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* tlb 3-10-00; remove updating of background model */
/* tlb 7-23-99; add nsites and nsites_obs to model; keep lambda <= 1 */
/* tlb 7-12-99; multiply pXijtheta by not_o_ij in e_step before calc of zij */
/* tlb 7-02-99; remove clobbering of theta by using logtheta[] in model */
/* tlb 6-28-99; add prior to lambda estimation using wnsites */
/* tlb 6-23-99; use same method for background counts for all models */
/* tlb 6-23-99; use average probability of strand/rc strand for background */
/* tlb 6-23-99; remove regularization of background */
/* tlb 6-22-99; precalculate background probability prior to inner loop */
/* tlb 6-22-99; combine zoops.c into this file */
/* tlb 6-22-99; change theta to obs in many places */
/* tlb 6-21-99; changed strand directions to compute 
    szik[j] = Pr(z_ij=1, s_ijk=1 | X_i, \theta)
    zi[j] = Pr(z_ij=1 | X_i, \theta) = sum_k=0^NDIR-1 szik[j]
*/
/* tlb 9-13-94; added different strand directions */

/**********************************************************************/
/*
	EM algorithm
*/
/**********************************************************************/

#include "meme.h"


#include "time.h"

/* added by M.H. */
/* functions for TCM M-step */
double NewtonRaphsonMethod(
  MODEL *model,			/* the model */
  DATASET *dataset		/* the dataset */
);

double computeFirstDerivative(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  double lambda			/* current lambda */
);

double computeSecondDerivative(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  double lambda         /* current lambda */
);



/**********************************************************************/
/*
  m_step 

	Do the M step of EM.
	Calculate the expected residue frequencies using the
	posterior probabilities of the site starts for each
	sequence.

	Time: O(n_samples*lseq*lseq)
*/
/**********************************************************************/
void m_step(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  PRIORS *priors,		/* the priors */
  double wnsites 		/* weight on prior on nsites */
)
{
  int i, j, k, ii, jj;
  THETA theta = model->theta;		/* theta of motif */
  THETA obs = model->obs;		/* observed frequencies of motif */
  int w = model->w; 			/* width of motif */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  int ndir = invcomp ? 2 : 1;		/* number of strands */
  int alength = dataset->alength;	/* length of alphabet */
  int n_samples = dataset->n_samples;	/* number of sequences in dataset */
  SAMPLE **samples = dataset->samples;	/* samples */
  double *back = dataset->back;		/* background frequencies */
  double q=0;				/* $Q$; sum of motif expectations */
  PTYPE ptype = priors->ptype;		/* type of prior */
  PriorLib *plib = priors->plib;	/* library for mixture priors */
  double total, total_obs;		/* scratch */
  double loglike0 = 0;			/* log-likelihood of sites under null */
  double loglike1;			/* log-likelihood of sites undr motif */
  double *p_count = priors->prior_count;/* pseudo counts */

  /* M step */

  /* initialize the observed frequency matrix to zero */
  for (i=0; i<w; i++) {
    for (j=0; j < alength; j++) obs(i, j) = 0; 
  }

  /* calculate the expected frequencies of residues in the different 
     site positions 
  */
  /* get the expected COUNTS, $c_k$ */
  for (i=0; i < n_samples; i++) { 		/* sequence */
    SAMPLE *s = samples[i];			/* array of samples */
    int lseq = s->length;			/* length this seq */
    double sw = s->sw;				/* sequence weight */
    double *lcb = s->logcumback;		/* cumulative backgnd. prob. */
    double qi = 0;				/* sum of z_ij */
    int last_j = lseq-w;			/* last site start */

    if (lseq < w) continue;			/* skip sequence; too short */

    /* get the counts on the two strands (DNA) */
    for (k=0; k<ndir; k++) {			/* strand direction */
      BOOLEAN ic = (k==1);			/* doing - strand */
      double *szik = s->sz[k];			/* sz_ik */
      
      for (j=0; j<=last_j; j++) {		/* site start */
        int off = ic ? lseq-w-j : j;		/* - strand offset from rgt. */
        char *res = ic ? s->resic+off : s->res+off;	/* integer sequence */
	double z = szik[j] * sw;		/* weight for motif */

	/* calculate: E[theta | X, z] */
        for (ii=0; ii<w; ii++) {		/* position in sequence */
	  int r = res[ii];
          if (r<alength) {			/* normal letter */
	    obs(ii, r) += z;			/* motif counts */
          } else { 				/* 'X': spread z over all */
            for (jj=0; jj<alength; jj++) obs(ii, jj) += z * back[jj];
          }
	} /* sequence position */

	/* 
	  calculate log-likelihood of sites under background model
	*/
        loglike0 += z * Log_back(lcb, j, w);

	qi += z;				/* sum of z_ij */
      } /* site start */
    } /* strand */

    q += qi;					/* $Q$; sum of q_i */
  }
  model->mll_0 = loglike0;			/* site log like. background */

  /* 
     M step for theta: convert COUNTS $c_k$ to FREQUENCIES $f_k$ 
     and use frequencies as estimate for theta--- $p_k^{(t+1)} = f_k$ 
  */
  /* convert counts to frequencies and regularize using pseudo-counts */
  for (i=0; i<w; i++) {				/* width */
    total = total_obs = 0.0;			/* total count for position i */
    /* get the total observed counts */
    for (j=0; j<alength; j++) total_obs += obs(i, j);
    /* get pseudo counts using prior */
    if (ptype == Dmix || ptype == Mega || ptype == MegaP)
      mixture_regularizer(obs[i], plib, p_count);
    /* adjust counts and total count by prior counts */
    for (j=0; j<alength; j++) {
      total += theta(i, j) = obs(i, j) + p_count[j];
    }
    /* normalize counts to probabilities */
    for (j=0; j<alength; j++) {
      obs(i, j) /= (total_obs ? total_obs : 1);	/* normalize to frequencies */
      theta(i, j) /= total;			/* normalize to probability */ 
    }
  } /* w */

  /* palindrome: enforce constraints */
  if (model->pal) {
    palindrome(theta, theta, w, alength);
    palindrome(obs, obs, w, alength);
  }

  /* compute log likelihood of sites under motif model */
  loglike1 = 0;
  for (i=0; i<w; i++) {
    for (j=0; j<alength; j++) {
      double f = obs(i, j);			/* letter frequency */
      if (f) loglike1 += f * LOG(f);		/* motif likelihood */
    } /* letter */
  } /* position in site */
  model->mll_1 = loglike1 *= q;			/* site log like. under motif */



  /* added by M.H. */
  /* finding the lambda for the next EM iteration is different for TCM model if secondary structure information is given */
  if ((dataset->secondaryStructureFilename != NULL) && (model->mtype == Tcm)) {

		double best_lambda = NewtonRaphsonMethod(model,dataset);		/* use Newton Raphson Method to find the best lambda */

		model->nsites_obs = q;
  		model->lambda_obs = MIN(1.0, best_lambda);
  		model->nsites = q*(1-wnsites) + model->psites*wnsites;

		/* lambda = (q(1-wnsites) + psites*wnsites) / n*m
			--> since lambda_obs = q/(n*m) it follows
			lambda = lambda_obs - lambda_obs*wnsites + (psites*wnsites/(n*m))
		*/
  		model->lambda = MIN(1.0, (model->lambda_obs - model->lambda_obs*wnsites + model->psites*wnsites/wps(dataset, w)) );
		
	}else{
  		/* M step for observed lambda */
		model->nsites_obs = q;
  		model->lambda_obs = MIN(model->nsites_obs / wps(dataset, w), 1.0);

  		/* M step for estimated lambda */
  		model->nsites = q*(1-wnsites) + model->psites*wnsites;
  		model->lambda = MIN(model->nsites / wps(dataset, w), 1.0);
	}
} /* m_step */



/* added by M.H. */
/**********************************************************************/
/*
  NewtonRaphsonMethod
		do Newton- Raphson method to find the root of the first derivative of term2 of the TCM M-step
		(finding the lambda that maximizes term2) 
  return lambda
*/
/**********************************************************************/
double NewtonRaphsonMethod(
  MODEL *model,			/* the model */
  DATASET *dataset		/* the dataset */
)
{
	int j;
	double lambda = ((model->lambda_obs > 1E-20) ? model->lambda_obs : model->lambda);		/* starting point (check if lambda_obs != 0) */	
   double min_lambda = 0; 										 	/* min lambda is 0 */
   double max_lambda = 1.0 / (double)(model->w - 1);	 	/* maximum lambda for a given w in case whole sequence consists of motifs */

	/* iterations */
	for (j=0; j<MAX_NEWTON_ITERATIONS; j++) {
		double fx1 = computeFirstDerivative(model, dataset, lambda);
		double fx2 = computeSecondDerivative(model, dataset, lambda);
		double lambda_new = lambda - (fx1/fx2);

		/* check if we run out of the boundaries
			if so use 0.5*(lambda - min_lambda) as the new lambda if we run into the negatives
			and use 0.5*(max_lambda - lambda) as the new lambda if lambda gets to big 
		*/
		if (lambda_new < min_lambda) {
			lambda_new = 0.5*(lambda - min_lambda);
		}else if (lambda_new > max_lambda) {
			lambda_new = 0.5*(max_lambda - lambda);
		}
		
		/* check if we have converged */
		if ( fabs(lambda_new - lambda) < 1E-10) {
			lambda = lambda_new;
			break;
		}
		lambda = lambda_new;
	}

	return lambda;
}	
		


/* added by M.H. */
/**********************************************************************/
/*
  computeFirstDerivative
		compute first derivative of term2 of the TCM M-step for a given lambda
  return value
*/
/**********************************************************************/
double computeFirstDerivative(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  double lambda
)
{
  int i, j;
  int w = model->w; 			/* width of motif */
  int n_samples = dataset->n_samples;	/* number of sequences in dataset */
  SAMPLE **samples = dataset->samples;	/* samples */

  double firstSum=0;				/* \sum \sum ((Z_ij-1)*m*sigma_ij) / (1- lambda*m*sigma_ij) */
  double sum_Z = 0;				/* \sum	\sum Z_ij */
  double val = 0;					/* return value */	

  for (i=0; i < n_samples; i++) { 	/* sequence */
    SAMPLE *s = samples[i];			/* array of samples */
    int lseq = s->length;				/* length this seq */
    double sw = s->sw;					/* sequence weight */
    int last_j = lseq-w;				/* last site start */
	 int m = last_j + 1;
	 
    if (lseq < w) continue;			/* skip sequence; too short */

    for (j=0; j<=last_j; j++) {		/* site start */
		double z = s->z[j] * sw;		/* weight for motif */

		double a = ((z-1) * m * s->sigma[j]);
		double b = (1 - (lambda * m * s->sigma[j]));
		firstSum += a / b;

		sum_Z += z; 	
	 }	
  }

  sum_Z /= lambda;			/* second sum is \sum\sum Z_ij / lambda */
  
  val = firstSum + sum_Z;	/* value of the first derivative */

  return val;	
}


/* added by M.H. */
/**********************************************************************/
/*
  computeSecondDerivative
		compute second derivative of term2 of the TCM M-step for a given lambda
  return value
*/
/**********************************************************************/
double computeSecondDerivative(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  double lambda
)
{
  int i, j;
  int w = model->w; 			/* width of motif */
  int n_samples = dataset->n_samples;	/* number of sequences in dataset */
  SAMPLE **samples = dataset->samples;	/* samples */

  double firstSum=0;				/* \sum \sum ((Z_ij-1)*m^2*sigma_ij^2) / ((1- lambda*m*sigma_ij)^2) */
  double sum_Z = 0;				/* \sum \sum Z_ij */
  double val = 0;					/* return value */	

  for (i=0; i < n_samples; i++) { 	/* sequence */
    SAMPLE *s = samples[i];			/* array of samples */
    int lseq = s->length;				/* length this seq */
    double sw = s->sw;					/* sequence weight */
    int last_j = lseq-w;				/* last site start */
	 int m = last_j + 1;
	 
    if (lseq < w) continue;			/* skip sequence; too short */

    for (j=0; j<=last_j; j++) {		/* site start */
		double z = s->z[j] * sw;		/* weight for motif */

		double a = ((z-1) * m*m * s->sigma[j]*s->sigma[j]);
		double b = (1 - (lambda * m * s->sigma[j]));
		b = b*b;
		firstSum += a / b;

		sum_Z += z; 	

	 }	
  }

  sum_Z /= (lambda*lambda);			/* second sum is \sum\sum Z_ij / lambda^2 */
  
  val = firstSum - sum_Z;				/* value of the second derivative */

  return val;	
}



/**********************************************************************/
/*
	e_step

	Do the E step of EM. 

	OOPS and ZOOPS models.

	Updates z array.

	Returns log Pr(X | theta).

	Time: O(n_samples*lseq*w)

	See notes 9/13/94
*/
/**********************************************************************/
double e_step(
  MODEL *model,                 /* the model */
  DATASET *dataset		/* the dataset */
)
{
  int i, j, k, ii;
  MTYPE mtype = model->mtype;		/* type of model */
  THETA logtheta1 = model->logtheta;	/* motif log(theta) */
  int w = model->w;			/* width of motif */
  int n_samples = dataset->n_samples;	/* number of sequences */
  BOOLEAN invcomp = model->invcomp;	/* use reverse complement strand, too */
  int ndir = invcomp ? 2 : 1;		/* number of strands */
  double log_sigma = log(1.0/ndir);	/* log \sigma */
  double lambda = (mtype==Zoops) ? model->lambda : 0;	/* lambda */
  double gamma = (mtype==Zoops) ? MIN(1.0,(lambda*wps(dataset,w))/n_samples) : 0;
  double log_1mgamma = (mtype==Zoops) ? LOG(1.0 - gamma) : 0;
  double logpX;				/* log Pr(X | \theta) */

  /* E step */

  convert_theta_to_log(model, dataset);	/* convert theta to log(theta) */

  /* calculate all the posterior offset probabilities */
  logpX = 0.0;				/* prob X given theta */
  for (i=0; i < n_samples; i++){	/* sequence */
    SAMPLE *s = dataset->samples[i]; 	/* sequence struct */
    int lseq = s->length;		/* length of the sequence */
    int m = lseq - w + 1;		/* number of possible sites */
    double *zi = s->z;			/* Pr(z_ij=1 | X_i, \theta) */ 
    double *not_o = s->not_o;		/* Pr(v_ij = 1) */
    double *lcb = s->logcumback;	/* cumulative background probability */
    double log_gamma =  (m && mtype==Oops) ? -LOG((double)m) : 0;/* log (1/m) */
    double log_lambda = (m && mtype==Zoops) ? LOG(gamma / m) : 0;/* exact */
    double log_pXitheta;		/* log Pr(X_i | theta) */

    if (lseq < w) continue;		/* skip sequence; too short */

    /* calculate probabilities for each strand direction */
    for (k=0; k<ndir; k++) {			/* strand */
      BOOLEAN ic = (k==1);			/* doing - strand */
      double *szik = s->sz[k];		/* Pr(X_i | z_ij=1, s_ijk=1, \theta) */

      /* calculate P(X_i | z_ij=1, \phi) */
      for (j=0; j<m; j++) {			/* site start */

			/* added by M.H. */
			/* use log sigma_ij instead of log 1/m if secondary structure information is given */
			/* NOTE: log_sigma is the prior for + or - strand --> here only + strand --> log_sigma = 0 */
	      if (dataset->secondaryStructureFilename != NULL) {
				if (mtype==Oops) {
		 	      log_gamma = s->log_sigma[j];
				}else{
		 	      log_lambda = LOG(gamma * s->sigma[j]);
				}
			}

        	/* sum_j=1^m sum_k=0^ndir-1 Pr(X_i | z_ij=1, s_ijk=1, \phi) = log Pr(X_i | z_ij=1, s_ijk=1 \phi) */
			double log_pXijtheta = LOG(not_o[j]) + ((mtype==Oops) ? log_sigma+log_gamma : log_sigma+log_lambda);

			int off = ic ? lseq-w-j : j; 	/* - strand offset from rgt. */
        	char *res = ic ? s->resic+off : s->res+off;	/* integer sequence */

			/* calculate the probability of positions outside of the site */
      	log_pXijtheta += lcb[lseq] - Log_back(lcb, j, w);

			/* calculate the probability of positions in the site */
			for (ii=0; ii<w; ii++) 
				log_pXijtheta += logtheta1(ii, (int) res[ii]);

        	/* set log szik to: 
          Pr(X_i | z_ij=1, s_ijk=1, \theta) \sigma (\gamma or \lambda)
        	*/
			szik[j] = log_pXijtheta;

        	/* set log z_ij to: 
          sum_k=0^ndir-1 Pr(X_i|z_ij=1,s_ijk=1,\theta)\sigma(\gamma or \lambda)
        	*/
        	zi[j] = !ic ? log_pXijtheta : LOG_SUM(zi[j], log_pXijtheta);
      } /* site start */  
    } /* strand */

    /* compute Pr(X_i | \phi) */
    log_pXitheta = (mtype==Oops) ? LITTLE : lcb[lseq] + log_1mgamma;
    for (j=0; j < m; j++) {			/* site start */
      log_pXitheta = LOG_SUM(log_pXitheta, zi[j]);
    }

    /* calculate logpX = \sum_{i=1}^n Pr(X_i | \phi) */
    logpX += log_pXitheta;

    /* sz_ijk and z_ij: normalize, delog and account for erasing
      1) Pr(z_ij=1, s_ijk=1 | X_i, \phi) \approx
           P(z_ij=1, s_ijk=1 | X_i, \phi) P(v_ij = 1) 
      2) P(z_ij=1 | X_i, \phi, V) = sum_k Pr(z_ij=1, s_ijk=1 | X_i, \phi)
    */
    for (k=0; k<ndir; k++) {		/* strand */
      double *szik = s->sz[k];		/* Pr(X_i | z_ij=1, s_ijk=1, \phi) */

      for (j=0; j<m; j++) {		/* site start */
        szik[j] = exp(szik[j] - log_pXitheta) * not_o[j];
        zi[j] = (k==0) ? szik[j] : zi[j] + szik[j];
      } /* site */

    } /* strand */

    for (j=m; j<lseq; j++) {		/* tail of sequence */
      zi[j] = 0;
    }

  } /* sequence */

  return logpX/log(2.0);
} /* e_step */

/**********************************************************************/
/*
	convert_theta_to_log

	Convert theta to log(theta) and include Pr(X).
*/
/**********************************************************************/
extern void convert_theta_to_log(
  MODEL *model,				/* the model */
  DATASET *dataset			/* the dataset */
)
{
  int i, j;
  THETA theta = model->theta;			/* theta */
  THETA logtheta = model->logtheta;		/* log(theta) */
  int w = model->w;				/* width */
  int alength = dataset->alength;		/* length of alphabet */
  double *back = dataset->back;			/* background freqs */

  for (i=0; i<w; i++) {				/* position */
    for (j=0; j<alength; j++) {			/* letter */
      logtheta(i, j) = LOGL(theta(i,j));
    }
    logtheta(i, j) = LOGL(back[j]);		/* Pr(X) */
  }
} /* convert_theta_to_log */

