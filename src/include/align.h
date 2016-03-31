/*
 * $Id: align.h,v 1.1.1.1 2005/07/29 18:34:59 nadya Exp $
 * 
 * $Log: align.h,v $
 * Revision 1.1.1.1  2005/07/29 18:34:59  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */


#ifndef ALIGN_H
#  define ALIGN_H 

#include <macros.h>
#include <user.h>
#include <logodds.h>
#include <hash_alph.h>

extern int align(
  int imotif,                           /* motif number */
  LOGODDS logodds,			/* log-odds array */
  int seqno,				/* sequence number (from 1) 
					   if <= 0, print alignment, 
					   otherwise, print .motif file */
  double threshold,			/* align sites above this */
  char *sample_name,			/* name of sample */
  char *eseq,				/* integer-coded sequence */
  BOOLEAN d[4],				/* strand directions to use */
  int lseq,				/* length of sequence */
  int w,				/* length of site */
  double *scores,			/* array to put scores in */
  FILE *outfile				/* stream for output */
);

#endif
