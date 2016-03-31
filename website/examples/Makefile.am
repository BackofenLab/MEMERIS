## $Id: Makefile.am,v 1.3 2005/09/15 21:35:28 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.3  2005/09/15 21:35:28  nadya
## # remove extra trailing slash
## #
## # Revision 1.2  2005/09/14 00:39:52  nadya
## # add mast sample sequences
## #
## # Revision 1.1.1.1  2005/07/25 23:21:55  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

if WEBSITE
exdir = $(WEBDIR)/examples
ex_DATA = $(EXAMPLES)
endif

EXAMPLES = \
	At.fa \
        sample-kabat.seq \
        sample-motif.dna 

EXTRA_DIST = $(EXAMPLES)
