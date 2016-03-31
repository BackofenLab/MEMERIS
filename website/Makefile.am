## $Id: Makefile.am,v 1.6 2005/10/13 20:11:42 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.6  2005/10/13 20:11:42  nadya
## # install Globals.pm with or without website. Need for perl scripts
## #
## # Revision 1.5  2005/10/04 23:44:26  nadya
## # make LOGS/unconditionally. Always need it for temp files that gets
## # created even for man pages.
## #
## # Revision 1.4  2005/08/25 00:29:40  nadya
## # fix typo
## #
## # Revision 1.3  2005/08/24 21:31:14  nadya
## # make LOGS world writeable
## #
## # Revision 1.2  2005/08/01 01:35:40  nadya
## # change location of LOGS to depend on prefix
## #
## # Revision 1.1.1.1  2005/07/31 20:30:51  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

if WEBSITE
SUBDIRS = html cgi-bin images scripts examples
else
SUBDIRS = scripts 
endif
