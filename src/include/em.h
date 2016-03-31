/*
 * $Id: em.h,v 1.1.1.1 2005/07/29 18:37:45 nadya Exp $
 * 
 * $Log: em.h,v $
 * Revision 1.1.1.1  2005/07/29 18:37:45  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef EM_H
#  define EM_H

/* 6-28-99 tlb; remove zoops and mcm e_steps */
/* 6-28-99 tlb; add wnsitesw */

extern void em(
  MODEL *model,			/* the model */
  DATASET *dataset 		/* the dataset */
);
extern void m_step(
  MODEL *model,			/* the model */
  DATASET *dataset,		/* the dataset */
  PRIORS *priors,		/* the priors */
  double wnsitesw                /* weight on prior on nsites */
);
extern double e_step(
  MODEL *model,                 		/* the model */
  DATASET *dataset              		/* the dataset */
);
extern double tcm_e_step(
  MODEL *model,					/* the model */
  DATASET *dataset				/* the dataset */
);
extern double like_e_step(
  MODEL *model,                 	 	/* the model */
  DATASET *dataset              		/* the dataset */
);
extern double discretize(
  MODEL *model,					/* the model */
  DATASET *dataset				/* the dataset */
);
extern void set_z (
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
);
extern void convert_theta_to_log(
  MODEL *model,					/* the model */
  DATASET *dataset 				/* the dataset */
);
#endif
