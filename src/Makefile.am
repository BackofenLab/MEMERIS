## $Id: Makefile.am,v 1.4 2005/12/15 00:33:58 nadya Exp $
##
## Process this file with automake to produce Makefile.in
##
## # $Log: Makefile.am,v $
## # Revision 1.4  2005/12/15 00:33:58  nadya
## # always install meme.bin and mast.bin
## #
## # Revision 1.3  2005/10/20 00:20:52  tbailey
## # *** empty log message ***
## #
## # Revision 1.2  2005/10/02 00:24:28  nadya
## # make adding libadd conditional in include switches . Now shows in
## # the compilation only when needed, noty by default
## #
## # Revision 1.1.1.1  2005/07/31 18:59:08  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AUTOMAKE_OPTIONS = no-dependencies


if WANT_STRLCPY
ADDDIR = libadd
ADDINC = -I./libadd
READSEQ_ADDSRC = $(ADDDIR)/strlcpy.c
endif

if SERVER 
SERVER_PROGS = meme-server mast-server 
if PARALLEL
SUBDIRS = . parallel
endif
endif

bin_PROGRAMS = alphtype getsize readseq \
               memeris mast.bin meme-client.bin \
               $(SERVER_PROGS) 

SEQ = hash.c hash_alph.c read_seq_file.c background.c
MOTIF = logodds.c motifs.c regress.c

alphtype_SOURCES = alphtype.c hash_alph.c
getsize_SOURCES = getsize.c $(SEQ)
getsize_LDADD = -lm
readseq_SOURCES = readseq.c ureadseq.c $(READSEQ_ADDSRC)
meme_server_SOURCES = meme-server.c readwrite.c
mast_server_SOURCES = mast-server.c readwrite.c
memeris_SOURCES = clock.c display.c dpalign.c em.c init.c mast-util.c \
               pssm-distr.c meme.c meme_util.c seq2theta.c starts.c prior.c \
               justlike.c \
               oops.c tcm.c banner.c gamma.c likelihood.c subseq7.c \
               discretize.c mp.c message.c $(MOTIF) $(SEQ) llr.c io.c logs.c
memeris_LDADD = -lm
mast_bin_SOURCES = mast.c clock.c diagram.c mast-util.c pssm-distr.c \
               $(MOTIF) $(SEQ) 
mast_bin_LDADD = -lm
meme_client_bin_SOURCES = meme-client.c readwrite.c

INCLUDES = -I./include $(ADDINC) 

EXTRA_DIST = include libadd parallel diagram.y 
