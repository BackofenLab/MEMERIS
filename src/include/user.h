/*
 * $Id: user.h,v 1.1.1.1 2005/07/29 19:12:25 nadya Exp $
 * 
 * $Log: user.h,v $
 * Revision 1.1.1.1  2005/07/29 19:12:25  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*
	User settable parameters
*/

#ifndef user_h
#define user_h

#define MSN	 24 		/* maximum length of sample name */
				/* MSN + 40 < PAGEWIDTH (see meme.h) */
#define MAXALPH  28		/* maximum length of alphabet + 1 for 'X' */
#define MAXG 101		/* maximum number of motifs + 1 */
#define MAXSITE 300		/* maximum length of a site */
#define MINSITES 2		/* minimum number of sites in valid motif */
#define LLR_RANGE 200		/* range of scaled LLR statistic */

#define MINCONS 0.2		/* Display 'X' as consensus if no letter f > */ 

/* minimum allowable motif width before shortening; 
   never make less than 2 or will crash! */
#define MIN_W 8
/* maximum allowable length before shortening */
#define MAX_W 50

#define MNAME 20		/* names of known motifs */
#define NMOTIFS MAXG		/* maximum number of known motifs */

#endif
