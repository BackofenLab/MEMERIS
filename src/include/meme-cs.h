/*
 * $Id: meme-cs.h,v 1.5.4.1 2006/01/31 18:22:56 nadya Exp $
 * 
 * $Log: meme-cs.h,v $
 * Revision 1.5.4.1  2006/01/31 18:22:56  nadya
 * rm OS-dependent declarations of
 * accept(), bind(), connect(), socket(), listen(), htons().
 * They are defined in system header files.
 *
 * Revision 1.5  2005/12/16 23:11:11  nadya
 * put MACOSX check back
 * add removed log entries
 *
 * Revision 1.4  2005/12/16 05:08:36  tbailey
 * Add support for CYGWIN operating system.
 *
 * Revision 1.3  2005/12/15 05:32:10  tbailey
 * getsize.c: Change getsize behavior: unless -nd given, duplicate
 * sequences are flagged  and NOT counted.
 * meme-cs.h: remove some Architecture checks to fix Cygwin version
 *
 * Revision 1.2  2005/10/19 20:33:57  nadya
 * add checking for MACOSX to prevent clash with OS-defined
 * accept(),bind(),conenct()
 *
 * Revision 1.1.1.1  2005/07/29 18:44:08  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/* client/server stuff for meme */

#ifndef MEME_CS_H
# define MEME_CS_H

#include "config.h"
#include "macros.h"
#include "user.h"

#undef MIN
#undef MAX

#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/time.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
/*#include <stdio.h>*/
#include <errno.h>
#include <fcntl.h>
#include "readwrite.h"

#ifdef __cplusplus
  #define CONST
  extern "C" {
#else
  #define CONST const
#endif

#ifdef __cplusplus
}
#endif

#endif
