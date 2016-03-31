## $Id: Makefile.am,v 1.15.4.2 2006/01/24 02:24:29 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.15.4.2  2006/01/24 02:24:29  nadya
## # EXT is obsoleted by using printf()
## #
## # Revision 1.15.4.1  2006/01/24 02:13:32  nadya
## # add handling of EXT variable
## #
## # Revision 1.15  2005/12/16 23:34:33  nadya
## # reorganize scripts separation for serer and client.
## #
## # Revision 1.14  2005/10/28 23:48:51  nadya
## # set MEME_LOGS and MEME_DB to user provided
## # need both environment settings when starting servers.
## #
## # Revision 1.13  2005/10/05 06:17:56  nadya
## # rm .txt.csh: rule and anything related.
## # It breaks making meme_config.csh on solaris if non-gnu make is used
## #
## # Revision 1.12  2005/10/02 05:56:49  nadya
## # rm ".csh" extension from scripts. Add ".bin" to the executables names.
## # add creation of meme_config.csh to source the environment in csh scripts
## #
## # Revision 1.11  2005/09/28 22:11:03  nadya
## # revert roc.csh to roc
## #
## # Revision 1.10  2005/09/23 20:07:06  nadya
## # change $< to $? in target rules. Some make variants on solaris
## # don't understand former.
## #
## # Revision 1.9  2005/09/21 20:37:42  nadya
## # add variable to specify user to run servers
## #
## # Revision 1.8  2005/08/31 04:51:55  nadya
## # create script and target for "make check"
## #
## # Revision 1.7  2005/08/25 20:18:02  nadya
## # use vavriable for setting mpi-lam direcotry
## #
## # Revision 1.6  2005/08/24 04:17:20  nadya
## # change convert2html from .pm to .pl
## # add machid handling
## #
## # Revision 1.5  2005/08/18 01:40:18  nadya
## # new runtests, initial revision
## #
## # Revision 1.4  2005/08/14 19:55:28  nadya
## # change $$@ target. Does not work with some versions of make.
## #
## # Revision 1.3  2005/08/11 18:53:09  nadya
## # rename from 2html. Perl does not understand modules if their name start as numeric.
## #
## # Revision 1.2  2005/08/09 23:08:30  nadya
## # simplify rules.
## #
## # Revision 1.1.1.1  2005/07/31 19:42:29  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

SEDSPEC = \
	-e 's%@WHICHPERL@%$(PERL)%' \
	-e 's%@PERLMODDIR@%$(perlmoddir)%' \
	-e 's%@MACH_ID@%$(MACHID)%' \
	-e 's%@sendmail@%$(SENDMAIL)%' \
	-e 's%@mail@%$(MAIL)%'\
	-e 's%@LAM_BIN@%$(LAMBIN)%'\
	-e 's%@MPICMD@%$(MPI_CMD)%' \
	-e 's%@MEMEDIR@%$(MEME_DIRECTORY)%' \
	-e 's%@MEMELOGS@%$(MEME_LOGS)%' \
	-e 's%@MEMEDB@%$(MEME_DB)%' \
	-e 's%@HOST_NAME@%$(HOSTNAME)%' \
	-e 's%@MASTPORT@%$(MAST_PORT)%'\
	-e 's%@MEMEPORT@%$(MEME_PORT)%' \
	-e 's%@MEMEUSER@%$(MEME_USER)%' \
	-e 's%@Q_TYPE@%$(QTYPE)%' \
	-e 's%@Q_NAME@%$(NAME)%' \
	-e 's%@MAX_TIME@%$(MAXTIME)%' \
	-e 's%@N_PROCS@%$(NPROCS)%' \
	-e 's%@MPI_TYPE@%$(MPITYPE)%'

CONFIG = meme_config meme_config.csh

if SERVER
SCRIPTS_SERVER = \
    lrun \
    mailer \
    memed \
    plot-usage \
    start-mast \
    start-meme \
    texit \
    update-plot-usage \
    usage-reports 
endif

SCRIPTS_CLIENT = \
    download \
    cat_max \
    compare_dates \
    fasta-get-markov \
    fasta-make-index \
    fasta-shuffle-letters \
    fasta-unique-names \
    make_logodds \
    mast \
    mast-client \
    mast2html \
    meme \
    meme-client \
    meme2html \
    plotgen \
    roc \
    rocket \
    runtests \
    sd \
    transfac2meme \
    update_mast_databases \
    update_meme_tests 

REST = \
    download \
    plot-usage \
    plotgen \
    rocket \
    sd \
    texit \
    usage-reports \
    update_meme_tests  

etcdir = $(prefix)/etc
etc_DATA = $(CONFIG)

perlmoddir=$(libdir)/perl
perlmod_SCRIPTS = convert2html.pl

bin_SCRIPTS = $(SCRIPTS_SERVER) $(SCRIPTS_CLIENT) 
TESTS = runcheck

convert2html.pl: convert2html.txt
	$(SED) $(SEDSPEC) $? > $@

meme2html: meme2html.pl
	$(SED) $(SEDSPEC) $? > $@

mast2html: mast2html.pl
	$(SED) $(SEDSPEC) $? > $@

TXTFILES = \
	cat_max.txt \
	compare_dates.txt \
	fasta-get-markov.txt \
	fasta-make-index.txt \
	fasta-shuffle-letters.txt \
	fasta-unique-names.txt \
	make_logodds.txt \
	lrun.txt \
	mailer.txt \
	mast-client.txt \
	meme-client.txt \
	meme.txt \
	mast.txt \
	meme_config.txt \
	meme_config.csh.txt \
	start-mast.txt  \
	start-meme.txt  \
	transfac2meme.txt \
	update-plot-usage.txt \
	update_mast_databases.txt \
	memed.txt \
   	runtests.txt \
   	roc.txt

TXT_FORMATTED = ${TXTFILES:.txt=}

cat_max: cat_max.txt 
	$(SED) $(SEDSPEC) $? > $@
compare_dates: compare_dates.txt 
	$(SED) $(SEDSPEC) $? > $@
fasta-get-markov: fasta-get-markov.txt 
	$(SED) $(SEDSPEC) $? > $@
fasta-make-index: fasta-make-index.txt 
	$(SED) $(SEDSPEC) $? > $@
fasta-shuffle-letters: fasta-shuffle-letters.txt 
	$(SED) $(SEDSPEC) $? > $@
fasta-unique-names: fasta-unique-names.txt 
	$(SED) $(SEDSPEC) $? > $@
lrun: lrun.txt 
	$(SED) $(SEDSPEC) $? > $@
make_logodds: make_logodds.txt 
	$(SED) $(SEDSPEC) $? > $@
mailer: mailer.txt 
	$(SED) $(SEDSPEC) $? > $@
mast-client: mast-client.txt 
	$(SED) $(SEDSPEC) $? > $@
meme-client: meme-client.txt 
	$(SED) $(SEDSPEC) $? > $@
mast: mast.txt 
	$(SED) $(SEDSPEC) $? > $@
meme: meme.txt 
	$(SED) $(SEDSPEC) $? > $@
start-mast: start-mast.txt  
	$(SED) $(SEDSPEC) $? > $@
start-meme: start-meme.txt  
	$(SED) $(SEDSPEC) $? > $@
transfac2meme: transfac2meme.txt
	$(SED) $(SEDSPEC) $? > $@
update-plot-usage: update-plot-usage.txt 
	$(SED) $(SEDSPEC) $? > $@
update_mast_databases: update_mast_databases.txt 
	$(SED) $(SEDSPEC) $? > $@
memed: memed.txt
	$(SED) $(SEDSPEC) $? > $@
runtests: runtests.txt
	$(SED) $(SEDSPEC) $? > $@
roc: roc.txt
	$(SED) $(SEDSPEC) $? > $@
meme_config: meme_config.txt
	$(SED) $(SEDSPEC) $? > $@
meme_config.csh: meme_config.csh.txt
	$(SED) $(SEDSPEC) $? > $@

CLEANFILES = \
	$(TXT_FORMATTED) \
	convert2html.pl \
	meme2html \
	mast2html 

EXTRA_DIST = \
	$(TXTFILES) \
	$(REST) \
	$(TESTS) \
	convert2html.txt \
	meme2html.pl \
	mast2html.pl 

