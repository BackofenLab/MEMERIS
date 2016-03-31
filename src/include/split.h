/*
 * $Id: split.h,v 1.1.1.1 2005/07/29 19:11:18 nadya Exp $
 * 
 * $Log: split.h,v $
 * Revision 1.1.1.1  2005/07/29 19:11:18  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef SPLIT_H
#  define SPLIT_H 

#include <macros.h>

/* longest line allowd in in file */
#define MAXLINE 10000

extern FILE *split(
  int nmode,		/* 0 - open; 1 - print */
  FILE *infile,
  int n
);

#endif
