/*
 * $Id: io.c,v 1.1.1.1.4.1 2006/01/26 09:16:27 tbailey Exp $
 * 
 * $Log: io.c,v $
 * Revision 1.1.1.1.4.1  2006/01/26 09:16:27  tbailey
 * Rename local function getline() to getline2() to avoid conflict with
 * function defined in stdio.h.
 *
 * Revision 1.1.1.1  2005/07/29 00:21:03  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

#include <macros.h>

/******************************************************************************/
/*
	getline2

	Read a newline or EOF-terminated line from a file.

	Returns a pointer to the line or NULL if at EOF when called.
*/
/******************************************************************************/
extern char *getline2(
  FILE *stream 				/* input stream */
) {
  char *s = NULL;			/* string to return */
  int c;
  int i = 0;				/* current position in string */
  while ((c=getc(stream)) != EOF) {
    if (i % GLBUFSIZ == 0) Resize(s, i+GLBUFSIZ, char); 
    s[i++] = c;
    if (c == '\n') break;		/* end of line */
  }
  if (feof(stream) && i==0) {
    return NULL;
  } else {
    return s;
  }
} /* getline2 */
