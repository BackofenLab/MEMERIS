/*
 * $Id: motifs.h,v 1.1.1.1 2005/07/29 18:45:20 nadya Exp $
 * 
 * $Log: motifs.h,v $
 * Revision 1.1.1.1  2005/07/29 18:45:20  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/* motifs.h */

#ifndef MOTIFS_H
#  define MOTIFS_H

extern int read_motifs (
  FILE *fdata,                          /* opened dataset file */
  char *filename,                       /* motif file */
  MOTIF motif[NMOTIFS],                 /* motif info */
  BOOLEAN save_dataset,                 /* return dataset in memory */
  DATASET *dataset                      /* integer-encoded dataset */
);

#endif
