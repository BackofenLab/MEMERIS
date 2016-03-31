## $Id: Makefile.am,v 1.9 2005/10/25 21:25:20 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.9  2005/10/25 21:25:20  nadya
## # make LOGS/ location configurable during "configure" step
## #
## # Revision 1.8  2005/10/13 20:11:42  nadya
## # install Globals.pm with or without website. Need for perl scripts
## #
## # Revision 1.7  2005/09/13 19:32:50  nadya
## # don't use reserved names
## #
## # Revision 1.6  2005/09/13 18:48:45  nadya
## # make clean target independent of conditional variable
## #
## # Revision 1.5  2005/08/20 02:21:49  nadya
## # add Validation module
## #
## # Revision 1.4  2005/08/11 17:41:29  nadya
## # rm trailing space
## #
## # Revision 1.3  2005/08/10 21:02:40  nadya
## # mv MAXTIME to Globals module
## #
## # Revision 1.2  2005/08/07 06:20:31  nadya
## # clean typos, rm debug
## #
## # Revision 1.1.1.1  2005/07/31 20:40:05  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = \
	-e 's%@SITE_CONTACT@%$(CONTACT)%'  \
	-e 's%@MEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@WEB_DIR@%$(WEBDIR)%' \
	-e 's%@MAX_TIME@%$(MAXTIME)%' \
	-e 's%@SITE_URL@%$(URL)%' 

SUFFIXES= .txt .pm
.txt.pm:
	$(SED) $(SEDSPEC) $< > $@

PM = Globals.pm Validation.pm
moduledir = $(libdir)/perl

if WEBSITE
module_SCRIPTS = $(PM)
else 
module_SCRIPTS = Globals.pm
endif

CLEANFILES = $(PM)
EXTRA_DIST = Globals.txt Validation.txt
