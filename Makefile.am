## $Id: Makefile.am,v 1.13 2005/12/16 23:29:08 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.13  2005/12/16 23:29:08  nadya
## # make install-data-hook target available in non-server mode
## #
## # Revision 1.12  2005/10/25 21:25:20  nadya
## # make LOGS/ location configurable during "configure" step
## #
## # Revision 1.11  2005/10/13 20:11:41  nadya
## # install Globals.pm with or without website. Need for perl scripts
## #
## # Revision 1.10  2005/10/13 19:58:53  nadya
## # add m4/ files implicitly
## #
## # Revision 1.9  2005/10/06 20:06:00  nadya
## # add definitions for configure macros that are missing from autoconf v < 2.54
## # will need to run bootstrap on install hostbefore running configure.
## #
## # Revision 1.8  2005/10/04 23:44:24  nadya
## # make LOGS/unconditionally. Always need it for temp files that gets
## # created even for man pages.
## #
## # Revision 1.7  2005/10/04 23:18:40  nadya
## # need etc/ always, contains init scripts.
## #
## # Revision 1.6  2005/09/16 00:22:05  nadya
## # put bootstrap back to top dir
## #
## # Revision 1.5  2005/09/15 21:29:24  nadya
## # add db directory creation
## #
## # Revision 1.4  2005/09/01 00:51:43  nadya
## # add test target for "make test"
## #
## # Revision 1.3  2005/08/24 21:57:32  nadya
## # clean autom4te.cache
## #
## # Revision 1.2  2005/08/01 00:25:15  nadya
## # change default tar name configuration, conflicts with rolls
## #
## # Revision 1.1.1.1  2005/07/25 21:33:48  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SUBDIRS = src scripts tests etc website

if SERVER
install-data-hook: dbdir
else
install-data-hook: 
endif

dbdir:
	mkdir -p $(MEME_DB)

install-data-local:
	mkdir -p $(MEME_LOGS)
	chmod a+w $(MEME_LOGS)

distdir = $(PACKAGE)_$(VERSION)

EXTRA_DIST = bootstrap m4/functions.m4 m4/headers.m4 Uniform_bfile GetSecondaryStructureValues.perl ReplaceUbyT.perl README_GetSecondaryStructureValues README_MEMERIS Changes_MH

distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	rm -rf autom4te.cache

test:
	cd scripts; $(MAKE) check
