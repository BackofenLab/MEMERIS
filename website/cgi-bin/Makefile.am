## $Id: Makefile.am,v 1.7 2005/09/13 19:32:50 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.7  2005/09/13 19:32:50  nadya
## # don't use reserved names
## #
## # Revision 1.6  2005/09/13 18:40:43  nadya
## # make clean target independent of conditional variable
## #
## # Revision 1.5  2005/08/19 22:59:20  nadya
## # move email validation to Validation module from check_address.cgi
## #
## # Revision 1.4  2005/08/19 01:25:44  nadya
## # use "mail" variable in meme.cgi and mast.cgi
## #
## # Revision 1.3  2005/08/10 21:02:39  nadya
## # mv MAXTIME to Globals module
## #
## # Revision 1.2  2005/08/10 20:33:24  nadya
## # use variable for MAXTIME, set by configure
## #
## # Revision 1.1.1.1  2005/07/25 23:17:00  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = \
        -e 's%@sendmail@%$(SENDMAIL)%' \
        -e 's%@mail@%$(MAIL)%'\
	-e 's%@WHICHPERL@%$(PERL)%'  \
	-e 's%@PERLLIBDIR@%$(libdir)/perl%' 


SUFFIXES= .pl .cgi
.pl.cgi:
	$(SED) $(SEDSPEC) $< > $@

CGI = job-status.cgi mast.cgi meme.cgi process_request.cgi

if WEBSITE
cgidir = $(WEBDIR)/cgi-bin
cgi_SCRIPTS = $(CGI)
endif

CLEANFILES = $(CGI)
EXTRA_DIST = job-status.pl mast.pl meme.pl process_request.pl
