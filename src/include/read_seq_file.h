/*
 * $Id: read_seq_file.h,v 1.1.1.1 2005/07/29 19:09:14 nadya Exp $
 * 
 * $Log: read_seq_file.h,v $
 * Revision 1.1.1.1  2005/07/29 19:09:14  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef READ_DATA_FILE_H
#define READ_DATA_FILE_H

extern DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  char *alpha,			/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  double seqfrac		/* fraction of input sequences to use */
);

extern BOOLEAN read_sequence(
  FILE *data_file,              /* file containing sequences */
  char **sample_name,           /* unique identifier of sequence */
  char **sample_de,             /* descriptor of sequence */
  char **sequence,              /* protein or DNA letters */
  long *length                   /* length of sequence */
);

#endif
