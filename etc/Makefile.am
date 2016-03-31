## $Id: Makefile.am,v 1.11 2005/12/13 00:58:10 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.11  2005/12/13 00:58:10  nadya
## # make creation of startup script conditional
## #
## # Revision 1.10  2005/10/25 21:25:20  nadya
## # make LOGS/ location configurable during "configure" step
## #
## # Revision 1.9  2005/10/13 21:15:51  nadya
## # move meme-explanation.html to etc/ to allow html conversion without installing web site
## #
## # Revision 1.8  2005/10/05 05:54:36  nadya
## # always put files in etc/. Needed for man-like output
## #
## # Revision 1.7  2005/09/26 23:15:16  nadya
## # add startup scripts for linux and solaris.
## #
## # Revision 1.6  2005/09/23 20:07:06  nadya
## # change $< to $? in target rules. Some make variants on solaris
## # don't understand former.
## #
## # Revision 1.5  2005/08/31 04:53:52  nadya
## # move usage info from mast.txt to mast.doc
## # add target mast.doc to Makefile.am
## #
## # Revision 1.4  2005/08/24 04:33:05  nadya
## # use MEME_DB for database location. Given as an option to configure
## #
## # Revision 1.3  2005/08/14 20:04:32  nadya
## # change $$@ target dependency. Does not work with some verisons of make
## #
## # Revision 1.2  2005/08/10 20:12:00  nadya
## # add environment setting scripts
## #
## # Revision 1.1.1.1  2005/07/29 23:37:56  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = \
	-e 's%@MEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' 

SHELLFILES_IN = meme.sh.in meme.csh.in
STARTUP_IN = meme.linux.in meme.solaris.in meme.irix.in
ETCFILES = \
    mast_db.list \
	mast_upstream_db.list \
    prior1.plib \
    prior30.plib \
    meme.doc \
    mast.doc \
    meme-explanation.html

if WANT_STARTUP
STARTUP = meme.$(MACHID)
endif

etcdir = $(prefix)/etc
etc_SCRIPTS = meme.sh meme.csh
etc_DATA = $(ETCFILES)

startupdir = $(prefix)/etc/startup
startup_SCRIPTS = $(STARTUP)

meme.sh: meme.sh.in
	$(SED) $(SEDSPEC) $? > $@

meme.csh: meme.csh.in
	$(SED) $(SEDSPEC) $? > $@

meme.$(MACHID): meme.$(MACHID).in
	$(SED) $(SEDSPEC) $? > $@

CLEANFILES = meme.sh meme.csh meme.$(MACHID)

EXTRA_DIST = $(SHELLFILES_IN) $(ETCFILES) $(STARTUP_IN)
