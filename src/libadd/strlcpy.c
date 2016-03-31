#include <stdio.h>
#include "config.h"


size_t strlcpy(char *dst, const char *src, size_t dstsize)
{
  int i;
  for (i=0; src[i] != '\0'; i++) {
    if (i<dstsize) dst[i] = src[i];
  }
  if (i<dstsize) dst[i] = '\0'; else dst[dstsize-1] = '\0';
  return(i); 
}

