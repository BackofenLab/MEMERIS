/*
 * $Id: io.h,v 1.1.1.1.4.1 2006/01/26 09:16:27 tbailey Exp $
 * 
 * $Log: io.h,v $
 * Revision 1.1.1.1.4.1  2006/01/26 09:16:27  tbailey
 * Rename local function getline() to getline2() to avoid conflict with
 * function defined in stdio.h.
 *
 * Revision 1.1.1.1  2005/07/29 18:40:41  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef IO_H
#define IO_H

extern char *getline2(
  FILE *stream 						/* input stream */
);
#endif
