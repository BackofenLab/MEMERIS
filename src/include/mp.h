/*
 * $Id: mp.h,v 1.1.1.1 2005/07/29 18:45:39 nadya Exp $
 * 
 * $Log: mp.h,v $
 * Revision 1.1.1.1  2005/07/29 18:45:39  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Bill Grundy 					       *
*								       *
***********************************************************************/
#ifndef MPMYID_H
#define MPMYID_H

extern void mpInit(
  int *argc,
  char **argv[]
);
extern void mpFinalize();
extern int mpNodes();
extern int mpMyID();
extern void mpReduceAdd(int *data);
extern void mpBroadcast(void *data, int bytes, int broadcaster);

#endif
