/*
 * $Id: starts.c,v 1.1.1.1.4.1 2006/01/24 21:08:38 nadya Exp $
 * 
 * $Log: starts.c,v $
 * Revision 1.1.1.1.4.1  2006/01/24 21:08:38  nadya
 * move to branch properly
 *
 * Revision 1.2  2006/01/20 03:03:32  tbailey
 * Remove obsolete (unused) argument sample_prob from call to subseq7().
 *
 * Revision 1.1.1.1  2005/07/29 17:27:03  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*                                                                      *
*       MEME++                                                         *
*       Copyright 1994, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                      *
*                                                                      *
***********************************************************************/

/* 7-29-99 tlb; change MAX_NSITES to be sites/w */
/* 6-29-99 tlb; add ic to subseq7 */
/*
	Routines to set up an array of starting points for EM.
*/

#include "meme.h"

static int get_nsites0 (
  int w,				/* width of motif */
  double min_nsites,			/* minimum nsites0 */ 
  double max_nsites,			/* maximum nsites0 */
  S_POINT *s_points			/* array of starting points */
);

/**********************************************************************/
/*
       	get_starts

	Create a list of starting points for EM. 

	If dataset->prob > 0, samples subsequences in the input
	dataset for good starting points.  
	Otherwise, e_cons is specified as the starting point.

	w is sampled in geometric progression by x w_factor;
	nsites is sampled in a geometric progression by x 2. 

	Returns number of starting points in list.

*/
/**********************************************************************/

extern S_POINT *get_starts(
  DATASET *dataset,		/* the dataset */
  MODEL *model,                 /* the model */
  char *e_cons,			/* encoded consensus sequence */
  double w_factor,		/* factor between sampled widths */
  int *n_starts			/* number of starting points */
)
{
  int i, j;
  int w;			/* width of motif */
  int n = 0;			/* number of starting points in array */
  S_POINT *s_points = NULL;	/* array of starting points */
  BOOLEAN pal = dataset->pal;	/* force DNA palindromes */
  THETA map = dataset->map; 	/* letter by frequency mapping matrix */
  double sample_prob=dataset->prob; /* sampling probability for subsequences */
  MTYPE mtype = model->mtype;	/* type of model */
  BOOLEAN ic = model->invcomp;	/* use reverse complement DNA strand, too */
  int min_w = model->min_w;	/* minimum width */
  int max_w = model->max_w;	/* maximum width */
  double min_nsites = model->min_nsites;	/* min. allowed sites */
  double max_nsites = model->max_nsites;	/* max. allowed sites */
  int incr, max_incr = pal ? 1 : 0;		/* do w, w+1 for palindromes */
  /* 
    try all values of w up to max_w in geometric progression 
  */
  for (w = min_w; w <= max_w; w = MIN(max_w, MAX(w+2, (int)(w*w_factor) ) ) ) {
  /*for (w = min_w; w <= max_w; w++) {	*/	/* all widths */

    for (incr=0; incr<=max_incr; incr++) {	/* w, w+1 for palindromes */
      int w0 = MIN(w+incr, MAXSITE);
      int n_nsites0 = max_nsites-min_nsites+1; 	/* upper bound on n_nsites0 */

      /* get list of nsites0 for this width and append to s_points list */
      Resize(s_points, n+n_nsites0, S_POINT);
      n_nsites0 = get_nsites0(w0, min_nsites, max_nsites, s_points+n);
      Resize(s_points, n+n_nsites0, S_POINT);

      /* fill in the starting points with the subsequence to use */
      if (!e_cons && sample_prob != 0) {   	/* sample subsequences */
	n_nsites0 = subseq7(mtype, ic, map, dataset, w0, n_nsites0, 
	  s_points+n);
      } else {					/* don't sample subsequences */
	for (i=n; i<n+n_nsites0; i++) {
	  s_points[i].e_cons0 = e_cons;		/* use the input consensus */
	  s_points[i].score = BIG;		/* tag as good starting point */
	}
      }

      /* set up human readable consensus sequence for starting point */
      for (i=n; i<n+n_nsites0; i++) {		/* consensus */
	char *e_cons0 = s_points[i].e_cons0;
	for (j=0; j<w0; j++)
	  s_points[i].cons0[j] = (e_cons0 ? unhash(e_cons0[j]) : 'x');
	s_points[i].cons0[j] = 0;
	if (PRINT_STARTS) {
	  printf("s=%d, score=%6.0f, w0=%3d, nsites0=%5.0f, cons=%s\n",
	    i, s_points[i].score, s_points[i].w0, 
	    s_points[i].nsites0, s_points[i].cons0);
	}
      } /* consensus */

      /* update length of starting point list */
      n += n_nsites0;

    } /* try w, w+1 for palindromes */

    /* end of w loop; tricky so max_w will be used */
    if (w == max_w) break;
  }

  *n_starts = n;			/* number of new starts */
  return s_points;
} /* get_starts */

/**********************************************************************/
/*	
	get_nsites0

	Get a list of the values to try for nsites0 and
	put them in the s_point array.   The array is set
	up with each entry with score LITTLE.

	List is geometric, increasing by factor of 2 from
	min_nsites to max_nsites.
	
	Returns the size of the added list.
*/
/**********************************************************************/
static int get_nsites0 (
  int w,				/* width of motif */
  double min_nsites,			/* minimum nsites0 */ 
  double max_nsites,			/* maximum nsites0 */
  S_POINT *s_points			/* array of starting points */
)
{
  double nsites;			/* number of sites */
  int n;				/* number of nsites for this w */

  /* initialize the starting points, making sure everything is initalized
     so that MPI won't barf
  */
  for (n=0, nsites=min_nsites; nsites < 2*max_nsites; n++, nsites*=2) {
  /*for (n=0, nsites=min_nsites; nsites < 2*max_nsites; n++, nsites+=1) {*/
    s_points[n].score = LITTLE;
    s_points[n].iseq = 0;
    s_points[n].ioff = 0;
    s_points[n].w0 = w;
    s_points[n].nsites0 = nsites<max_nsites ? nsites : max_nsites;
    s_points[n].wgt_nsites = 0;
    s_points[n].e_cons0 = NULL;
    s_points[n].cons0 = NULL; 
    Resize(s_points[n].cons0, w+1, char);
    s_points[n].cons0[0] = '\0';
    if (nsites >= max_nsites) {n++; break;}
  }

  return n;				/* number of starts for w */
} /* get_nsites0 */
