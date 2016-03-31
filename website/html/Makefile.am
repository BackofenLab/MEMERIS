## $Id: Makefile.am,v 1.7.4.1 2006/01/24 21:14:57 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.7.4.1  2006/01/24 21:14:57  nadya
## # create COPYRIGHT.htmlk from new header file and copyright notice
## #
## # Revision 1.7  2005/12/21 20:20:07  nadya
## # add meme-install.html
## #
## # Revision 1.6  2005/10/13 21:15:52  nadya
## # move meme-explanation.html to etc/ to allow html conversion without installing web site
## #
## # Revision 1.5  2005/09/01 22:14:34  nadya
## # create template files to use version automatically
## #
## # Revision 1.4  2005/08/25 00:48:28  nadya
## # check for existing file before installing
## #
## # Revision 1.3  2005/08/24 19:22:17  nadya
## # add hook to make a symbolic link for index.html
## #
## # Revision 1.2  2005/08/24 05:29:06  nadya
## # add job-status and general-faq files
## #
## # Revision 1.1.1.1  2005/07/31 20:22:38  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = \
	-e 's%@version@%$(VERSION)%' 

TMPL = head-copyright
INFILES = intro.in meme.in mast.in 
IN_FORMATTED = ${INFILES:.in=.html}

SUFFIXES = .in .html
.in.html:
	$(SED) $(SEDSPEC) $< > $@

COPYRIGHT.html: $(TMPL) $(top_srcdir)/COPYING
	/bin/cat $(TMPL) $(top_srcdir)/COPYING > COPYRIGHT.html
	echo "</PRE>" >> COPYRIGHT.html
	echo "</BODY>" >> COPYRIGHT.html
	echo "</HTML>" >> COPYRIGHT.html

if WEBSITE
htmldir = $(WEBDIR)
html_DATA = $(FILES) $(IN_FORMATTED) COPYRIGHT.html
install-data-hook:
	cd $(WEBDIR); \
	if [ ! -e index.html ]; then \
	    ln -s intro.html index.html; \
	fi
endif

FILES = \
	general-faq.html \
	help_alphabet.html \
	help_distribution.html \
	help_format.html \
	job-status.html \
	mast-databases.html \
	mast-faq.html \
	mast-input.html \
	mast-intro.html \
	mast-output-example.html \
	mast-output.html \
	meme-download.html \
	meme-faq.html \
	meme-input.html \
	meme-intro.html \
	meme-install.html \
	meme-output-example.html \
	motif-format.html \
	papers.html \
	release_notes.html 

CLEANFILES = $(IN_FORMATTED) COPYRIGHT.html

EXTRA_DIST = $(FILES) $(INFILES) $(TMPL)

