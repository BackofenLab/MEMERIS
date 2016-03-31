## $Id: Makefile.am,v 1.5 2005/10/25 21:22:09 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.5  2005/10/25 21:22:09  nadya
## # add pssm-distr.c to list of sources
## #
## # Revision 1.4  2005/10/07 02:48:45  nadya
## # create a rule for making an *.o file with mpicc
## # use $< in place of $? to compile on SUN
## #
## # Revision 1.3  2005/10/07 02:43:14  nadya
## # *** empty log message ***
## #
## # Revision 1.2  2005/10/01 23:55:05  nadya
## # redefine make's macros COMPILE and LINK in order to compile with
## # correct MPI flags when using mpicc
## #
## # Revision 1.1.1.1  2005/07/31 19:23:53  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies

# redefine make macros 
CC = $(MPICC)
AM_CFLAGS = $(MPIFLAGS)
INCLUDES = -I../include $(MPIINC)
COMPILE = $(CC) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(AM_CFLAGS)
LINK = $(CCLD) $(AM_CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@

bin_PROGRAMS = meme_p
meme_p_SOURCES = \
        ../background.c \
        ../banner.c \
        ../clock.c \
        ../display.c \
        ../discretize.c \
        ../dpalign.c \
        ../em.c  \
        ../gamma.c \
        ../hash.c \
        ../hash_alph.c \
        ../init.c  \
        ../io.c \
        ../justlike.c \
        ../likelihood.c \
        ../llr.c \
        ../logs.c \
        ../logodds.c \
        ../mast-util.c \
        ../meme.c \
        ../meme_util.c \
        ../message.c \
        ../motifs.c \
        ../mp.c \
        ../oops.c \
        ../prior.c  \
        ../pssm-distr.c  \
        ../read_seq_file.c \
        ../regress.c \
        ../seq2theta.c \
        ../starts.c \
        ../subseq7.c \
        ../tcm.c 

# need this rule for making object files with mpicc
COMPILE_RULE = $(COMPILE) -c `test -f $? || echo '$(srcdir)/'`$?

background.o: ../background.c 
	$(COMPILE_RULE)
banner.o: ../banner.c 
	$(COMPILE_RULE)
clock.o: ../clock.c
	$(COMPILE_RULE)
display.o: ../display.c
	$(COMPILE_RULE)
discretize.o: ../discretize.c 
	$(COMPILE_RULE)
dpalign.o: ../dpalign.c
	$(COMPILE_RULE)
em.o: ../em.c
	$(COMPILE_RULE)
gamma.o: ../gamma.c 
	$(COMPILE_RULE)
hash.o: ../hash.c 
	$(COMPILE_RULE)
hash_alph.o: ../hash_alph.c 
	$(COMPILE_RULE)
init.o: ../init.c
	$(COMPILE_RULE)
io.o: ../io.c 
	$(COMPILE_RULE)
justlike.o: ../justlike.c 
	$(COMPILE_RULE)
likelihood.o: ../likelihood.c 
	$(COMPILE_RULE)
llr.o: ../llr.c 
	$(COMPILE_RULE)
logs.o: ../logs.c 
	$(COMPILE_RULE)
logodds.o: ../logodds.c 
	$(COMPILE_RULE)
mast-util.o: ../mast-util.c
	$(COMPILE_RULE)
meme.o: ../meme.c 
	$(COMPILE_RULE)
meme_util.o: ../meme_util.c 
	$(COMPILE_RULE)
message.o: ../message.c 
	$(COMPILE_RULE)
motifs.o: ../motifs.c 
	$(COMPILE_RULE)
mp.o: ../mp.c 
	$(COMPILE_RULE)
oops.o: ../oops.c 
	$(COMPILE_RULE)
prior.o: ../prior.c  
	$(COMPILE_RULE)
pssm-distr.o: ../pssm-distr.c  
	$(COMPILE_RULE)
read_seq_file.o: ../read_seq_file.c 
	$(COMPILE_RULE)
regress.o: ../regress.c
	$(COMPILE_RULE)
seq2theta.o: ../seq2theta.c 
	$(COMPILE_RULE)
starts.o: ../starts.c 
	$(COMPILE_RULE)
subseq7.o: ../subseq7.c 
	$(COMPILE_RULE)
tcm.o: ../tcm.c 
	$(COMPILE_RULE)

