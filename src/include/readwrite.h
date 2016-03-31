/*
 * $Id: readwrite.h,v 1.1.1.1 2005/07/29 19:09:33 nadya Exp $
 * 
 * $Log: readwrite.h,v $
 * Revision 1.1.1.1  2005/07/29 19:09:33  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#ifndef READWRITE_H
 #define READWRITE_H

#include <unistd.h>

extern int readline(register int fd,
		    register char *ptr,
		    register int maxlen);

extern int readn (register int fd,
		  register char *ptr,
		  register int nbytes);

extern int writen(register int fd,
		  register char *ptr,
		  register int nbytes);

#endif

