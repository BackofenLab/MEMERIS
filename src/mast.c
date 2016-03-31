/*
 * $Id: mast.c,v 1.4.4.2 2006/01/31 20:51:38 nadya Exp $
 * 
 * $Log: mast.c,v $
 * Revision 1.4.4.2  2006/01/31 20:51:38  nadya
 * rm initialization for 'name'
 *
 * Revision 1.4.4.1  2006/01/31 20:22:48  nadya
 * init name with zeros
 *
 * Revision 1.4  2005/10/25 19:02:26  nadya
 * change c++ style comment to proper c
 *
 * Revision 1.3  2005/10/20 00:20:52  tbailey
 * *** empty log message ***
 *
 * Revision 1.2  2005/10/01 23:58:04  nadya
 * update documentation comment
 *
 * Revision 1.1.1.1  2005/07/29 18:06:43  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/************************************************************************
*                                                                     	*
*       MAST                                                           	*
*       Author: Timothy L. Bailey                                      	*
*                                                                       *
*	Copyright							*
*	(1994 - 2001) The Regents of the University of California.	*
*	All Rights Reserved.						*
*									*
*	Permission to use, copy, modify, and distribute any part of 	*
*	this software for educational, research and non-profit purposes,*
*	without fee, and without a written agreement is hereby granted, *
*	provided that the above copyright notice, this paragraph and 	*
*	the following three paragraphs appear in all copies.		*
*									*
*	Those desiring to incorporate this software into commercial 	*
*	products or use for commercial purposes should contact the 	*
*	Technology Transfer Office, University of California, San Diego,*
*	9500 Gilman Drive, La Jolla, California, 92093-0910, 		*
*	Ph: (619) 534 5815.						*
*									*
*	IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO 	*
*	ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR 	*
*	CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF 	*
*	THE USE OF THIS SOFTWARE, EVEN IF THE UNIVERSITY OF CALIFORNIA 	*
*	HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 		*
*									*
*	THE SOFTWARE PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE *
*	UNIVERSITY OF CALIFORNIA HAS NO OBLIGATIONS TO PROVIDE 		*
*	MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.  *
*	THE UNIVERSITY OF CALIFORNIA MAKES NO REPRESENTATIONS AND 	*
*	EXTENDS NO WARRANTIES OF ANY KIND, EITHER EXPRESSED OR IMPLIED, *
*	INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 	*
*	MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE, OR THAT 	*
*	THE USE OF THE MATERIAL WILL NOT INFRINGE ANY PATENT, 		*
*	TRADEMARK OR OTHER RIGHTS.  					*
************************************************************************/

/**********************************************************************/
/*
  mast <logodds> <database> <alphabet> [optional arguments...]
	See <install-path>/bin/mast for documentation.
	See <install-path>/bin/make_logodds for format of logodds file.
*/
/**********************************************************************/

#define DEFINE_GLOBALS
#include "mast.h"
#include "diagram.h"
#include "projrel.h"

/* printing */
BOOLEAN debug = FALSE;		/* turn on debugging features */

#ifndef EXP
#define EXP 0
#else
#define EXP 1
#endif

/* size of memory chunks for sortable sequences */
#define RCHUNK 100

/* default e-thresh */
#define EXPECT 10

/* maximum pairwise motif correlation */
#define MAXCORR 0.6

/* constants for printing control */
#define MMSN 34
#define PAGEWIDTH 80
#define MAXID (80 - MMSN - 7 - PVALUE - LENGTH)
#define PVALUE 8
#define LENGTH 6
#define MAXIDLINES 10

/* number of different integral score values */
#define RANGE 100

/* macro to copy one open file to another; noop if file NULL */
#define copy_file(in_file, out_file) \
  if (in_file != NULL) { \
    int c; \
    rewind(in_file); \
    while ((c = getc(in_file)) != EOF) putc(c, out_file); \
  }

/*
  Data Definitions
*/

/* sortable sequence score record */
typedef struct {
  long  fp;		/* file pointer to beginning of sequence rec. */
  long length;		/* length of sequence */
  double *score;	/* motif scores */
  int *loc;		/* motif locations */
  int pos;		/* positive motif 0 score */
  int neg;		/* negative motif 0 score */
  double Pvalue;	/* p-value of product of p-values */
  BOOLEAN known_pos;	/* known positive if true (for ROC) */
  char *id;		/* name of sequence */
  int strand;		/* -1 neg. strand, 0 both/protein, +1 pos. strand */
  double *comp;		/* actual sequence composition */
  double **pv;		/* pvalue distribution for each motif */
} SORT_SEQ;

/*
  Local Procedures
*/

static int get_motifs(
  BOOLEAN xlate_dna,		/* database is DNA and motifs protein */
  char *nlogodds,		/* name of log-odds matrix */
  char *alphabet,		/* alphabet of logodds matrices */
  char *blast_alphabet,		/* corresponding BLAST alphabet */
  int *p[MAXASCII],		/* alphabet permutation/substitution matrix */
  BOOLEAN shuffle,		/* shuffle the motif columns */
  double min_ev,		/* minimum E-value of motifs to use */
  int umotifs, 			/* number motifs to be used */
  int mlist[MAXLO],		/* list of motifs given in -m option */
  BOOLEAN motifs[MAXLO],	/* list of motifs to output */
  LO *los[MAXLO],		/* array of logodds matrices */
  double *f,			/* array of null letter probabilities */
  int range 			/* set logodds matrices to [0..range] */
);

static SORT_SEQ *get_scores(
  BOOLEAN dna,		/* database is DNA */
  STYPE stype,		/* handling of DNA strands */
  BOOLEAN xlate_dna,	/* database is DNA and motifs protein */
  double f[],		/* null letter probability distribution */
  FILE *fdata,		/* the database */
  FILE *fsave,		/* good sequences if database on stdin */
  LO *los[],		/* array of pointers to log-odds matrices */
  int nmotifs,		/* number of log-odds matrices in los */
  int range, 		/* set entries in matrices to [0..range] */
  double **pv,		/* p-value tables for each motif */
  BOOLEAN sonly,	/* calculate p-value of observed spacings only */
  BOOLEAN lump,		/* combine spacings into one p-value */
  BOOLEAN use_seq_comp,	/* adjust E-value for actual sequence composition */
  double e_max,		/* maximum sequence E-value to print */
  int kmotifs,		/* number of known motifs */
  MOTIF motif[],	/* data returned by read_motifs */
  BOOLEAN status,	/* show progress */
  int min_seqs,		/* lower bound on nseqs */
  int *nseqs,		/* total sequences in database */
  double *residues,	/* total number of residues in seqs */
  int *n_hits		/* sequences less than e_max */
);

static double calc_roc(
  SORT_SEQ *seq,			/* array of score values */
  int nseqs,				/* number of sequences */
  int pos,				/* number of known positives */
  int roc_num,				/* number of fp's to truncate ROC at */
  double thresh				/* threshold for error calculation */
);

static int p_compare(
  const void *v1,
  const void *v2
);

static double calc_p_value(
  char *sample_name,	/* name of sample */
  long length,           /* length of the sequence */
  STYPE stype,		/* handling of DNA strands */
  int nmotifs,          /* number of motifs */
  LO *los[],            /* array of pointers to log-odds matrices */
  double *best_score,   /* array of best score for each motif */
  int *best_loc,        /* position of best match for each motif */
  double range,		/* number of different score values */
  double **pv,		/* p-value tables for each motif */
  BOOLEAN sonly,	/* calculate p-value of observed spacings only */
  BOOLEAN lump,		/* combine spacings into one p-value */
  int norder,		/* number of motifs in diag */
  BOOLEAN debug 	/* turn on debugging features */
);

static void print_mast_results(
  BOOLEAN read_stdin,			/* read database on stdin */
  char *alphabet,			/* alphabet */
  double *back,				/* letter probability distribution */
  char *bfile,				/* no background file */
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of different strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN use_seq_comp,			/* adjust E-value for sequence composition */
  FILE* hit_file,			/* table of hits */
  FILE* diag_file,			/* table of motif diagrams */
  FILE* note_file,			/* table of annotated sequences */
  BOOLEAN doc,				/* print documentation */
  int rank,				/* rank of first result to print */
  double e_max,				/* maximum sequence E-value to print */
  double m_thresh,			/* maximum motif p-value to print */
  double w_thresh,	 		/* maximum motif p-value--weak hits */
  BOOLEAN use_seq_p,			/* use sequence not position p-values */
  int n_hits,				/* sequences less than e_max */
  LO *los[],				/* array of pointers to lo matrices */
  int bad_list[], 			/* list of highly correlated motifs */
  int nbad,				/* number of bad motifs */
  BOOLEAN rem_corr,			/* remove correlated motifs */
  char *mfname, 			/* motif file name to print */
  int nmotifs,				/* number motifs read */ 
  char *database,			/* name of database */
  char *dfname,				/* name of database to print */
  int nseqs, 				/* number of sequences in database */
  double res				/* total number of residues in seqs */
);

static int make_mast_tables(
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of DNA strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN use_seq_comp,			/* adjust E-value for sequence composition */
  BOOLEAN best_motifs, 			/* diagrams have only best motif */
  int kmotifs,	 			/* number of known motifs */
  MOTIF motif[],			/* data returned by read_motifs */
  char *motif_file,			/* motif file name */
  FILE *fdata,				/* the database file */
  int nseqs,				/* number of sequences scored */
  SORT_SEQ *seqs,			/* sortable array of seq/score records*/
  int n_hits,				/* size of seqs array */
  FILE* hit_file,			/* table of hits */
  FILE* diag_file,			/* table of motif diagrams */
  FILE* note_file,			/* table of annotated sequences */
  int rank,				/* rank of first result to print */
  int smax,				/* maximum sequences to print */
  double e_max,				/* maximum E-value to print */
  double m_thresh,			/* maximum motif p-value to print */
  double w_thresh,	 		/* maximum motif p-value-- weak hits */
  BOOLEAN use_seq_p,			/* use sequence not position p-values */
  BOOLEAN hit_list,                     /* create hit list instead of diagram */
  LO *los[],				/* array of pointers to lo matrices */
  int nmotifs 				/* number motifs read */ 
);

static void print_annotation(
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of DNA strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  double thresh,			/* threshold for strong hits */
  FILE *note_file,			/* table of annotated sequences */
  LO *los[],				/* array of pointers to lo matrices */
  char *sequence,			/* sequence of sample */
  long length,				/* length of sample */
  char *sample_name,			/* name of sample */
  char *id,				/* identifier text for sample */
  int strand,				/* strand of DNA or 0 if both/protein */
  double pvalue,			/* combined p-value of sequence */
  double evalue,			/* combined E-value of sequence */
  TILING tiling				/* tiling and diagram of sequence */
);

static char *best_frame(
  int nmotifs, 				/* number of motifs */
  long length,				/* length of sequence */
  TILING				/* tiling of sequence */
);

static int get_bad_motifs(
  double maxcorr,			/* maximum allowed corellation */
  int nmotifs,				/* number of motifs */ 
  LO *los[],				/* array of logodds matrices */
  int bad_list[] 			/* list of highly correlated motifs */
);

/**********************************************************************/
/*
  main routine
*/
/**********************************************************************/
extern int main (int argc, char *argv[])
{
  /* defaults */
  char *nlogodds = "";
  char *database = "";
  char *alphabet = "";
  char *mfname = NULL;			/* don't print the motif file name */
  char *dfname = NULL;			/* print actual name of database */
  BOOLEAN print_seq = TRUE;		/* print annotated sequence */
  BOOLEAN weak = FALSE;			/* print weak hits */
  int rank = 1;				/* print results starting with rank */
  int smax = 0;				/* print all results */
  double e_max = EXPECT;		/* maximum sequence p-value to print */
  double m_thresh = 1e-4;		/* maximum motif p-value to print */
  double w_thresh;			/* maximum motif p-value--weak hits 
					   default: m_thresh * 10 */
  double min_ev = BIG;			/* minimum E-value of motifs to use */
  int min_seqs = 0;			/* lower bound on number db sequences */
  int roc_num = 0;			/* don't compute ROC */
  BOOLEAN use_seq_p = FALSE;		/* use sequence p-value for m_thresh */
  BOOLEAN shuffle = FALSE;		/* don't shuffle the motif columns */
  BOOLEAN sonly = FALSE;		/* print product of spacing p-values */
  BOOLEAN lump = FALSE;			/* combine spacings into one p-value */
  BOOLEAN status = TRUE;		/* show progress */
  BOOLEAN doc = TRUE;			/* print documentation */
  STYPE stype = Combine;		/* handling of DNA strands */
  BOOLEAN xlate_dna = FALSE;		/* don't tranlate DNA sequences */
  BOOLEAN best_motifs = FALSE;          /* only print best motif in diagrams */
  BOOLEAN use_seq_comp = FALSE;		/* adjust E-values for actual compos. */
  BOOLEAN rem_corr = FALSE;		/* remove overly corellated motifs */ 
  BOOLEAN read_stdin = FALSE;		/* read database from stdin */ 
  BOOLEAN hit_list = FALSE;		/* print block diagrams */ 
  
  /* locals */
  int i, j, k;
  char *arguments;			/* saved argv[1] */
  FILE *fdata = NULL;			/* the database */
  FILE *fsave = NULL;			/* good seqs if reading stdin */
  FILE *hit_file, *diag_file, *note_file;/* temporary files for tables */
  char *motif_file = NULL;		/* motif file name */
  char *bfile = NULL;			/* no background file */
  int nseqs = 0;			/* total sequences in database */
  double residues = 0;			/* total number of residues in seqs */
  int n_hits = 0;			/* sequences less than e_max */
  BOOLEAN dna;				/* true if sequences are DNAB */
  BOOLEAN error;			/* true if something wrong */

  /* motifs and logodds matrices */
  int nmotifs;				/* number motifs read */ 
  int umotifs=0;			/* number motifs to be used */
  int kmotifs=0;	 		/* number of known motifs */
  LO *los[MAXLO];			/* array of logodds matrices */
  BOOLEAN motifs[MAXLO];		/* list of motifs to output */
  int bad_list[MAXLO]; 			/* list of highly correlated motifs */
  int nbad;				/* number of bad (correlated) motifs */
  int mlist[MAXLO];			/* list of motifs given in -m option */
  MOTIF motif[NMOTIFS];			/* data returned by read_motifs */
  DATASET dummy;			/* dummy arg for read_motifs */
  char *blast_alphabet;			/* corresponding BLAST alphabet */
  int *perm[MAXASCII];			/* permutation/substitution matrix */

  /* things for normalization */
  SORT_SEQ *seqs;		/* sortable array of seq/score records */
  double **pv = NULL;		/* p-value tables for each motif */
  double *back;			/* background letter probability distribution */

#ifdef MALLOC_DEBUG
  int malloc_debug(int);
  malloc_debug(2);
#endif

  (void) myclock();             	/* record CPU time */

  /* save script arguments, set command to blank and shift arguments */
  arguments = argv[1]; argv[1] = ""; argv++; argc--;

  /* get the command line arguments */
  i = 1;
  DO_STANDARD_COMMAND_LINE(1,
    NON_SWITCH(1, \r,
      switch (i++) {
        case 1: break;	/* dummy */
        default: COMMAND_LINE_ERROR;
      }
    );
    DATA_OPTN(2, logodds, , , nlogodds = _OPTION_);
    DATA_OPTN(2, database, , , database = _OPTION_);
    DATA_OPTN(2, alphabet, , , alphabet = _OPTION_);
    FLAG_OPTN(2, stdin, , read_stdin = TRUE);
    FLAG_OPTN(1, sep,
      \t\tscore reverse complement DNA strand as a separate \n\t\t\tsequence, 
      stype = Separate);
    FLAG_OPTN(1, norc, \t\tdo not score reverse complement DNA strand, 
      stype = Norc);
    FLAG_OPTN(1, dna, \t\ttranslate DNA sequences to protein, xlate_dna = TRUE);
    FLAG_OPTN(1, comp, \t\tadjust p-values and E-values for sequence composition, 
      use_seq_comp = TRUE);
    DATA_OPTN(1, rank, <rank>,
      \tprint results starting with <rank> best (default: 1),
      rank = atoi(_OPTION_));
    DATA_OPTN(1, smax, <smax>,
      \tprint results for no more than <smax> sequences\n\t\t\t(default: all),
      smax = atoi(_OPTION_));
    DATA_OPTN(1, ev, <ev>,
      \tprint results for sequences with E-value < <ev>\n\t\t\t(default: 10),
      e_max = atof(_OPTION_)); 
    DATA_OPTN(1, mt, <mt>, 
      \tshow motif matches with p-value < mt (default: 0.0001), 
      m_thresh = atof(_OPTION_)); 
    FLAG_OPTN(1, w, 
      \t\tshow weak matches (mt<p-value<mt*10) in angle brackets, weak = TRUE);
    DATA_OPTN(1, bfile, <bfile>, \tread background frequencies from <bfile>,
      bfile = _OPTION_);
    FLAG_OPTN(1, seqp, \t\tuse SEQUENCE p-values for motif thresholds\n\t\t\t(default: use POSITION p-values), use_seq_p = TRUE);
    DATA_OPTN(1, mf, <mf>, \tprint <mf> as motif file name, mfname = _OPTION_);
    DATA_OPTN(1, df, <df>, \tprint <df> as database name, dfname = _OPTION_);
    DATA_OPTN(1, minseqs, <minseqs>, \tlower bound on number of sequences in db,
      min_seqs = atoi(_OPTION_));
    DATA_OPTN(1, mev, <mev>,+\tuse only motifs with E-values less than <mev>, 
      min_ev = atof(_OPTION_));
    DATA_OPTN(1, m, <m>,+\tuse only motif(s) number <m> (overrides -mev), 
      mlist[umotifs++] = atoi(_OPTION_));
    DATA_OPTN(1, diag, <diag>, \tnominal order and spacing of motifs,
      diagram = _OPTION_);
    FLAG_OPTN(1, best, \t\tinclude only the best motif in diagrams,
      best_motifs = TRUE);
    FLAG_OPTN(1, remcorr, \tremove highly correlated motifs from query,
      rem_corr = TRUE);
    FLAG_OPTN(1, brief, 
      \tbrief output--do not print documentation, doc = FALSE);
    FLAG_OPTN(1, b, \t\tprint only sections I and II, print_seq = FALSE);
    FLAG_OPTN(1, nostatus, \tdo not print progress report, status = FALSE);
    FLAG_OPTN(1, hit_list, \tprint hit_list instead of diagram; implies -text, hit_list = TRUE);
    /* DEBUG OPTIONS */
    FLAG_OPTN(EXP, deb, \t\tprint p-values and exit; implies -text, debug = TRUE);
    DATA_OPTN(EXP, k, <k>, 
      \tfile of known motifs to mark as "+" or for ROC if \n\t\t\t-r given, 
      motif_file = _OPTION_);
    DATA_OPTN(EXP, r, <r>, \tcompute roc_<r> using <k>.tag file, 
      roc_num = atoi(_OPTION_));
    FLAG_OPTN(EXP, shuffle, \tshuffle columns of motifs, shuffle = TRUE);
    FLAG_OPTN(EXP, sonly, \tuse only spacing p-values in product, sonly = TRUE);
    FLAG_OPTN(EXP, lump, \t\tcombine spacings into one p-value, lump = TRUE);
  );

  /* check input */
  if (rank < 1) {
    fprintf(stderr, "<rank> must be at least 1.\n");
    exit(1); 
  }

  /* set e_max to BIG if doing ROC or printing all sequence p-values */
  if (roc_num > 0 || debug) e_max = BIG;

  /* set maximum p-value for weak hits */
  if (weak) {
    w_thresh = m_thresh * 10.0;			/* 10 times less likely */
  } else {
    w_thresh = m_thresh;			/* no weak threshold */
  }

  /* open the datafile or standard input */
  if (!read_stdin) {
    if (!(fdata = fopen(database, "r"))) {
      fprintf(stderr, "Cannot open file `%s'.\n", database);
      exit(1); 
    }
  } else {					/* reading stdin */
    fdata = stdin;
    fsave = tmpfile();				/* save good sequences here */
  } /* open database or standard input */

  /* initialize the background frequencies and alphabets */
  back = init_mast_background(bfile, alphabet, stype, xlate_dna, 
    &blast_alphabet, perm);

  /* get the motifs and parse the ordering and spacing diagram */
  nmotifs = get_motifs(xlate_dna, nlogodds, alphabet, blast_alphabet, perm,
    shuffle, min_ev, umotifs, mlist, motifs, los, back, RANGE);

  /* get list of motifs that should be removed */
  nbad = get_bad_motifs(MAXCORR, nmotifs, los, bad_list);

  /* remove highly correlated motifs from query */
  if (rem_corr && nbad) {
    for (i=j=k=0; i<nmotifs; i++) {
      if (los[i]->name < bad_list[k] || k >= nbad) {
        los[j++] = los[i];
      } else {
        k++;
      }
    }
    nmotifs -= nbad;
  }

  /*
    Set current sequence alphabet to DNA if translating.
  */
  if (xlate_dna) setalph(0);

  /* determine the type of database being searched */
  dna = (xlate_dna || los[0]->dna);

  /* get the type of handling of DNA strands */
  if (!dna) {					/* protein database */
    if (stype == Separate || stype == Norc) {
      fprintf(stderr, 
        "You may not specify -sep or -norc with a protein database.\n");
      exit(1);
    }
    stype = Protein;
  } /* database type */

  /* read in (optional) known motif occurrences */
  if (motif_file) kmotifs = read_motifs(fdata, motif_file, motif, 0, &dummy);

  /*
    calculate p-values of all integer score values in range [0...w*RANGE] 
  */
  Resize(pv, nmotifs, double *);
  for (i=0, error=FALSE; i<nmotifs; i++) {
    /*pv[i] = calc_cdf(los[i], RANGE, back);*/
    pv[i] = calc_pssm_cdf(los[i]->w, los[i]->alen, RANGE, los[i]->logodds, back);
    if (!pv[i]) {
      fprintf(stderr, "There is something wrong with motif %d\n", los[i]->name);
      error = TRUE;
    }
  }
  if (error) exit(1);			/* error with motifs */

  /* 
     get the scores and p-values for all sequences in the database and
     sort them by p-value in ascending order 
  */
  seqs = get_scores(dna, stype, xlate_dna, back, fdata, fsave, los, nmotifs, 
    RANGE, pv, sonly, lump, use_seq_comp, e_max,
    kmotifs, motif, status, min_seqs, &nseqs, &residues, &n_hits
  );

  /* print p-values and exit if debug on */
  if (debug) {
    for (i=1; i<nmotifs; i++) {                 /* from motif */
      int j;
      printf("CORR  %2d ", los[i]->name);
      for (j=0; j<i; j++) {                     /* to motif */
        printf(" %5.2f", los[i]->corr[j]);
      } /* to motif */
      printf("\n");
    } /* from motif */
    for (i=0; i<n_hits; i++) {
      printf("NORM %g\n", seqs[i].Pvalue);
    } /* hit */
    exit(0);
  }

  /* calculate the ROC of the motifs and exit if doing ROC */
  if (roc_num > 0) {
    printf("\n\nROC %d = %f\n", roc_num, 
      calc_roc(seqs, n_hits, motif[0].pos, roc_num, los[0]->thresh));
    exit(0);
  }

  /* create the three output tables in temporary files */
  hit_file = tmpfile();
  diag_file = tmpfile();
  note_file = print_seq ? tmpfile() : NULL;
  if (fsave) fdata = fsave;			/* saved sequences */
  n_hits = make_mast_tables(dna, stype, xlate_dna, use_seq_comp, best_motifs, 
    kmotifs, motif, motif_file, fdata, nseqs, seqs, n_hits, hit_file, 
    diag_file, note_file, rank, smax, e_max, m_thresh, w_thresh, use_seq_p, 
    hit_list, los, nmotifs
  );

  /* print the results */
  print_mast_results(
    read_stdin, alphabet, back, bfile, dna, stype, xlate_dna, use_seq_comp,
    hit_file, diag_file, note_file, doc, 
    rank, e_max, m_thresh, w_thresh, use_seq_p,
    n_hits, los, bad_list, nbad, rem_corr, mfname, nmotifs, database, dfname,
    nseqs, residues
  );

  /* print the command line to the script which was saved in argv[1] */
  printf("%s\n", arguments);

  /* finish status report */
  if (status) fprintf(stderr, "\n");

  return 0;
} /* main */

/**********************************************************************/
/*
        calc_roc
 
        Calculate the ROC (receiver operating characteristic) of
        a sorted list of score values with known positives flagged.
*/
/**********************************************************************/
static double calc_roc(
  SORT_SEQ *seqs,			/* array of score values */
  int nseqs,                            /* number of sequences */
  int pos,				/* number of known positives */
  int roc_num,				/* number of fp's to truncate ROC at */
  double thresh				/* threshold for error calculation */
)
{
  int i;
  double roc = 0;			/* receiver operating characteristic */
  double tpp = 0; 			/* true positive proportions */
  double fpp = 0;			/* false positive proportions */
  double tp = 0;			/* true positives so far */
  double fp = 0;			/* false positives so far */
  double newtpp, newfpp;
  int neg = nseqs - pos;		/* number of known negatives */
  printf("thresh = %f\n", thresh);

  /* loop over sequences, best score first */
  for (i=0; i<nseqs; i++) {		/* sequence */

    /* update tp and fp since score has changed (or last score reached) */
    /* known positive? */
    if (seqs[i].known_pos) {
      tp++;
    } else {
      fp++;
    }

    /* trapezoidal rule : (y2 - y1)/2 dx */
    newtpp = MIN(1.0, tp / pos);	/* seqs marked "?" in prosite are
					   marked but not counted in pos
					   so newtpp could go over 1.0 */
    newfpp = fp / neg;
    roc += .5 * (newtpp + tpp) * (newfpp - fpp); 
    tpp = newtpp;
    fpp = newfpp;
 
    if (seqs[i].known_pos) { printf("+ "); } else { printf("- "); }
    printf("%-*.*s %8.1e\n", MMSN, MMSN, seqs[i].id, seqs[i].Pvalue);

    /* reached limiting number of false positives? */
    if (fp >= roc_num) break;

  } /* sequence */

  /* normalize by fpp to get ROC <roc_num> */
  if (fpp == 0) {
    roc = 1.0;
  } else {
    roc /= fpp;
  }

  return roc;
} /* calc_roc */

/**********************************************************************/
/*
        p_compare
 
        Compare two p-values in ascending order.  Return <0 >0
        if the second is >, < the first.  If they are equal,
	resolves ties by returning <0 if the second has smaller fp.
*/
/**********************************************************************/
static int p_compare(
  const void *v1,
  const void *v2
)
{
  const SORT_SEQ * s1 = (const SORT_SEQ *) v1;
  const SORT_SEQ * s2 = (const SORT_SEQ *) v2;
  double diff = s2->Pvalue - s1->Pvalue;
  if (diff == 0) diff = (double) (s1->fp - s2->fp);
  return ((diff > 0) ? -1 : ( (diff < 0) ? 1 : 0) );
} /* p_compare */

/**********************************************************************/
/*
	calc_p_value

   	Calculate the p-value of the product of the individual motif 
	p-values < observed

	Returns the combined p-value.
*/
/**********************************************************************/
static double calc_p_value(
  char *sample_name,	/* name of sample */
  long length,           /* length of the sequence */
  STYPE stype,		/* handling of DNA strands */
  int nmotifs,          /* number of motifs */
  LO *los[],            /* array of pointers to log-odds matrices */
  double *best_score,   /* array of best score for each motif */
  int *best_loc,        /* position of best match for each motif */
  double range,		/* number of different score values */
  double **pv,		/* p-value tables for each motif */
  BOOLEAN sonly,	/* calculate p-value of observed spacings only */
  BOOLEAN lump,		/* combine spacings into one p-value */
  int norder,		/* number of motifs in diag */
  BOOLEAN debug 	/* turn on debugging features */
)
{
  int i;
  double pvalue;
  int nspaces;			/* number of motif pairs spacings given for */
  double k_scores;		/* product of motif score p-values */
  double k_spacing;		/* product of spacing p-values */
  int smotifs;			/* number of motifs that got scored */

  /* calculate the product of the motif score p-values */
  k_scores = 1.0;
  for (i=smotifs=0; i<nmotifs; i++) {
    int ws = los[i]->ws;	/* width of motif in sequence */
    int x = best_score[i];	/* observed EV (extreme value) */
    long n = length - ws + 1;	/* number of samples in EV */
    double p;			/* p-value of EV */
    if (best_score[i] == LITTLE) continue;	/* motif wasn't scored */
    if (stype == Combine) n*=2;	/* combining both DNA strands */
    EV(pv[i][x], n, p);		/* compute sequence p-value of x */
    k_scores *= p;		/* product of p-values */
    smotifs++;			/* number of motifs scored */
  }

  /* multiply by p-values of motif spacings if provided (norder > 1) */
  k_spacing = 1.0;
  nspaces = 0;			/* number of spaces p-valued */
  for (i=1; i<norder; i++) {	/* space to previous motif */
    if (space[i] >= 0) {	/* don't ignore this space */
      double p;			/* to hold p-value of obs. spacing */
      int err;			/* error in pos ith and i-1th motifs */
      int mi = order[i];	/* current motif */
      int mim1 = order[i-1];	/* previous motif */ 
      long npos;		/* number of positions for motif */
      int ws_avg;		/* average width of adjacent motifs */
      long mind;		/* minimum spacing that fits seq */
      long maxd;		/* maximum spacing that fits seq */

      /* get absolute error: diffence between observed and nominal spacing */
      err = abs(space[i] - (best_loc[mi] - best_loc[mim1]));

      /* get the average motif width and maximum number of placements */
      ws_avg = (int) (los[mim1]->ws + los[mi]->ws)/2.0;	/* round down */
      npos = length - ws_avg + 1;

      mind = MAX(space[i] - err, ws_avg - length);
      maxd = MIN(space[i] + err, length - ws_avg);
      if (mind >= 0) {
	p = (maxd-mind+1) * (npos - (maxd+mind)/2.0);
      } else {
	p = -mind * (npos - (1-mind)/2.0) + (maxd+1) * (npos - maxd/2.0);
      }
      p /= npos * npos;
      if (p > 1.0 || p < 0.0) fprintf(stderr, 
	"\nerror in spacing p-value:%8.8s L=%8ld mind=%8ld E=%4d %-10g\n", 
	 sample_name, length, mind, err, p);

      /* skip if no possible spacing */
      if (maxd < mind) {
	continue;
      }

      k_spacing *= p;		/* product of p-values */
      nspaces++;		/* number of spacings used */

    } /* don't ignore */
  } /* space to previous motif */

  /* finish the calculation */
  if (sonly) {					/* spacings only */
    pvalue = qfast(nspaces, k_spacing);
  } else if (lump && nspaces) {			/* lump spacings */
    pvalue = qfast(nspaces, k_spacing);		/* spacing pvalue */
    pvalue = qfast(smotifs+1, k_scores*pvalue); /* spacing and motif p-value */
  } else {					/* spacings and motif scores */
    pvalue = qfast(smotifs+nspaces, k_scores*k_spacing);
  }

  return pvalue;
} /* calc_p_value */

/**********************************************************************/
/*
	print_mast_results

	Print MAST results.
*/
/**********************************************************************/
static void print_mast_results(
  BOOLEAN read_stdin,			/* read database on stdin */
  char *alphabet,			/* alphabet */
  double *back,				/* letter probability distribution */
  char *bfile,				/* no background file */
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of different strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN use_seq_comp,			/* compensate for seq composition */
  FILE* hit_file,			/* table of hits */
  FILE* diag_file,			/* table of motif diagrams */
  FILE* note_file,			/* table of annotated sequences */
  BOOLEAN doc,				/* print documentation */
  int rank,				/* rank of first result to print */
  double e_max, 			/* maximum sequence E-value to print */
  double m_thresh,			/* maximum motif p-value to print */
  double w_thresh,	 		/* maximum motif p-value--weak hits */
  BOOLEAN use_seq_p,			/* use sequence not position p-values */
  int n_hits,				/* sequences less than e_max */
  LO *los[],				/* array of pointers to lo matrices */
  int bad_list[], 			/* list of highly correlated motifs */
  int nbad,				/* number of bad motifs */
  BOOLEAN rem_corr,			/* remove correlated motifs */
  char *mfname, 			/* motif file name */
  int nmotifs,				/* number motifs read */ 
  char *database,			/* name of database */
  char *dfname,				/* name of database to print */
  int nseqs, 				/* number of sequences in database */
  double res 				/* total number of residues in seqs */
)
{
  int i, j;
  char *database_date;			/* creation/modif. date of database */
  struct stat stbuf;			/* buffer for stat call */
  char *stars = 
"********************************************************************************";
  char *mtype = (xlate_dna || !dna) ? "peptide" : "nucleotide";
  char *dbtype = (!dna) ? "peptide" : "nucleotide";

  /* 
    announce the program 
  */
  i = strlen(ARCHIVE_DATE) - 9;
  printf("%s\n", stars);
  printf("MAST - Motif Alignment and Search Tool\n");
  printf("%s\n", stars);
  printf(
"\tMAST version %s (Release date: %*.*s)\n\n"
"\tFor further information on how to interpret these results or to get\n"
"\ta copy of the MAST software please access http://meme.nbcr.net.\n",
    VERSION, i, i, ARCHIVE_DATE+7);
  printf("%s\n", stars);

  /* 
    print reference citation 
  */
  printf("\n\n%s\n", stars);
  printf("REFERENCE\n");
  printf("%s\n", stars);
  printf(
"\tIf you use this program in your research, please cite:\n"
"\n"
"\tTimothy L. Bailey and Michael Gribskov,\n"
"\t\"Combining evidence using p-values: application to sequence homology\n"
"\tsearches\", Bioinformatics, 14(48-54), 1998.\n"
  );
  printf("%s\n", stars);

  /* 
    print info on the database 
  */
  printf("\n\n%s\n", stars);
  printf("DATABASE AND MOTIFS\n");
  printf("%s\n", stars);
  printf("\tDATABASE %s (%s)\n", dfname ? dfname : database, dbtype);
  /* print name (and date if running in UNIX) of database */
#ifdef UNIX
  if (!read_stdin) {		/* Get date of database if not reading stdin */
    stat(database, &stbuf);
    database_date = ctime(&stbuf.st_mtime);
    printf("\tLast updated on %s", database_date);
  }
#endif
  /* print number of sequences in database */
  i = stype==Separate ? 2 : 1;			/* only count + strand */
  printf("\tDatabase contains %d sequences, %.0f residues\n\n", nseqs/i, res/i);
  /* print handling of DNA strands */
  if (stype == Combine) {
    printf("\tScores for positive and reverse complement strands are combined.\n\n");
  } else if (stype == Separate) {
    printf("\tPositive and reverse complement strands are scored separately.\n\n");
  } else if (stype == Norc) {
    printf("\tReverse complement strands are not scored.\n\n");
  }


  /* 
    print info on the motifs 
  */
  if (mfname != NULL) printf("\tMOTIFS %s (%s)\n", mfname, mtype);
  printf("\tMOTIF WIDTH BEST POSSIBLE MATCH\n");
  printf("\t----- ----- -------------------\n");
  for (i=0; i<nmotifs; i++) {
    printf("\t%3d   %3d   %*.*s\n", los[i]->name, los[i]->w, 
      los[i]->w, los[i]->w, los[i]->best_seq);
  }
  if (diagram != NULL) 
    printf("\n\tNominal ordering and spacing of motifs:\n\t %s\n", diagram);

  /* 
   print motif correlations 
  */
  if (nmotifs > 1) {
    printf("\n\tPAIRWISE MOTIF CORRELATIONS:\n");
    printf("\tMOTIF");
    for (i=0; i<nmotifs-1; i++) printf(" %5d", los[i]->name);
    printf("\n");
    printf("\t-----");
    for (i=0; i<nmotifs-1; i++) printf(" -----");
    printf("\n");
    for (i=1; i<nmotifs; i++) {                 /* from motif */
      printf("\t  %2d ", los[i]->name);
      for (j=0; j<i; j++) {                     /* to motif */
        printf(" %5.2f", los[i]->corr[j]);
      } /* to motif */
      printf("\n");
    } /* from motif */
    if (nbad) {
      printf(
"\tCorrelations above %4.2f may cause some combined p-values and\n"
"\tE-values to be underestimates.\n", MAXCORR
      );
      if (rem_corr) {
        printf("\tRemoved motif");
      } else {
        printf("\tRemoving motif");
      }
      if (nbad > 1) printf("s "); else printf(" ");
      for (i=0; i<nbad; i++) {
        if (i==nbad-1) {
          if (i>0) printf(" and ");
        } else {
          if (i>0) printf(", ");
        }
        printf("%d", bad_list[i]);
      }
      if (rem_corr) {
        printf(
          " because they have correlation > %4.2f\n\twith the remaining motifs.\n", 
          MAXCORR);
      } else {
        printf(" from the query may be advisable.\n");
      }
    } else {
      printf("\tNo overly similar pairs (correlation > %4.2f) found.\n",
        MAXCORR);
    }
  } /* nmotifs > 1 */

  /* print background model frequencies */
  if (!use_seq_comp) {
    int i, pcol;
    char *c;
    printf("\n\tRandom model letter frequencies (from %s):",
      bfile ? bfile : "non-redundant database");
    for (c=alphabet, i=0, pcol=80; *c; c++, i++) {
      int lpos = xlate_dna ? protbhash(alphabet[i]) : hash(alphabet[i]);
      pcol += 8;          			/* start of printed thing */
      if (pcol >= 80) {pcol=15; printf("\n\t");}
      printf("%c %5.3f ", *c, back[lpos]);
    }
    printf("\n");
  } else {
    printf("\n\tUsing random model based on each target sequence composition.\n");
  }
 
  /* end database and motif section */
  printf("%s\n", stars);

  /* 
    print table of hits documentation
  */
  printf("\n\n%s\nSECTION I: HIGH-SCORING SEQUENCES\n%s\n", stars, stars);
  if (doc) {
    char *fs1 = xlate_dna && stype==Separate ? 
      "\n\t- The strand and frame of the (best) motif match(es) is shown." : 
      xlate_dna ? "\n\t- The frame of the (best) motif match(es) is shown." : "";
    char *fs2 = xlate_dna ? 
      "\n\t  Frames 1, 2, and 3 are labeled a, b c, respectively." : "";
    printf(
"\t- Each of the following %d sequences has E-value less than %g.\n",
      n_hits, e_max);
    if (rank > 1) {
      printf(
"\t- The %d best-matching sequences have been omitted.\n", rank-1);
    }
    printf(
"\t- The E-value of a sequence is the expected number of sequences\n"
"\t  in a random database of the same size that would match the motifs as\n"
"\t  well as the sequence does and is equal to the combined p-value of the\n"
"\t  sequence times the number of sequences in the database.\n"
"\t- The combined p-value of a sequence measures the strength of the\n"
"\t  match of the sequence to all the motifs and is calculated by\n"
"\t    o finding the score of the single best match of each motif\n"
"\t      to the sequence (best matches may overlap),\n");
    printf(
"\t    o calculating the sequence p-value of each score,\n"
"\t    o forming the product of the p-values,\n");
    if (norder > 1) {
      printf(
"\t    o multiplying by the p-value of the observed spacing of\n"
"\t      pairs of adjacent motifs (given the nominal spacing),\n");
    }
    printf(
"\t    o taking the p-value of the product.\n");
    printf(
"\t- The sequence p-value of a score is defined as the\n"
"\t  probability of a random sequence of the same length containing\n"
"\t  some match with as good or better a score.\n");
    printf(
"\t- The score for the match of a position in a sequence to a motif\n"
"\t  is computed by by summing the appropriate entry from each column of\n"
"\t  the position-dependent scoring matrix that represents the motif.\n"
"\t- Sequences shorter than one or more of the motifs are skipped.%s%s\n"
"\t- The table is sorted by increasing E-value.\n", fs1, fs2);
    printf("%s\n\n", stars);
  } /* doc */

  /* 
    print table of hits headings and data
  */
  {
    int il = xlate_dna || stype==Separate ? MAXID-3 : MAXID+4;/* length of id */
    char *st1 = xlate_dna ? "FRAME  " : stype==Separate ? "STRAND " : "";
    char *st2 = xlate_dna ? "-----  " : stype==Separate ? "------ " : "";
    char *f = "%-*s %-*s%s %8s %6s\n";
    printf(f, MMSN, "SEQUENCE NAME", il, "DESCRIPTION", st1,"E-VALUE ","LENGTH");
    printf(f, MMSN, "-------------", il, "-----------", st2,"--------","------");
    copy_file(hit_file, stdout);
    printf("\n%s\n\n", stars);
  }

  /* 
    print table of diagrams documentation
  */
  printf("\n\n%s\nSECTION II: MOTIF DIAGRAMS\n%s\n", stars, stars);
  if (doc) {
    char *fs0 = stype==Separate ? "the strand and\n\t  " : "";
    char *fs1 = stype==Combine ? "s" : "";
    char *fs2 = xlate_dna ? "f" : "";			
    char *fs3 = xlate_dna ? 
      "\n\t\t    in frame f.  Frames 1, 2, and 3 are labeled a, b c." : ".";
    char *fs4 = stype==Combine ? "\t\t    A minus sign indicates that the occurrence is on the\n\t\t    reverse complement strand.\n" : "";
    char *ptype = use_seq_p ? "SEQUENCE" : "POSITION";
    char *set = use_seq_p ?
      "some random subsequence in a set of n,\n where n is the sequence length minus the motif width plus 1," :
      "a single random subsequence of the length of the motif";
    printf(
"\t- The ordering and spacing of all non-overlapping motif occurrences\n"
"\t  are shown for each high-scoring sequence listed in Section I.\n"
"\t- A motif occurrence is defined as a position in the sequence whose\n"
"\t  match to the motif has %s p-value less than %g.\n"
"\t- The %s p-value of a match is the probability of\n"
"\t  %s\n"
"\t  scoring at least as well as the observed match.\n"
"\t- For each sequence, all motif occurrences are shown unless there\n"
"\t  are overlaps.  In that case, a motif occurrence is shown only if its\n"
"\t  p-value is less than the product of the p-values of the other\n"
"\t  (lower-numbered) motif occurrences that it overlaps.\n"
"\t- The table also shows %sthe E-value of each sequence.\n"
"\t- Spacers and motif occurences are indicated by\n"
"\t   o -d-    `d' residues separate the end of the preceding motif \n"
"\t\t    occurrence and the start of the following motif occurrence\n",
      ptype, w_thresh, ptype, set, fs0);
    printf(
"\t   o [%sn%s]  occurrence of motif `n' with p-value less than %g%s\n%s",
      fs1, fs2, m_thresh, fs3, fs4);
    if (w_thresh != m_thresh) printf(
"\t   o <%sn%s>  occurrence of motif `n' with %g < p-value < %g%s\n%s",
      fs1, fs2, m_thresh, w_thresh, fs3, fs4);
    printf("%s\n\n", stars);
  } /* doc */

  /*
    print table of diagrams headings and data
  */
  {
    char *st1 = stype==Separate ? "STRAND" : "";	/* strand heading */
    char *st2 = stype==Separate ? "------" : "";	/* strand underline */
    int nl = (*st1=='\0') ? MMSN : MMSN-4;		/* length of name */
    char *f = "%-*s%s %8s  %s\n";			/* format */
    printf(f, nl, "SEQUENCE NAME", st1, "E-VALUE ", "MOTIF DIAGRAM");
    printf(f, nl, "-------------", st2, "--------", "-------------");
    copy_file(diag_file, stdout);
    printf("\n%s\n\n", stars);
  }

  /* 
   print table of annotated sequences documentation and data
  */
  if (note_file) {
    printf("\n\n%s\n", stars);
    printf("SECTION III: ANNOTATED SEQUENCES\n");
    printf("%s\n", stars);
    if (doc) {
      char *fs1 = xlate_dna ? "/frame" : "";		/* frame string 3 */
      char *fs2 = stype==Combine ? 
          " (a minus sign indicates that\n\t  the occurrence is on the reverse complement strand)" : 
          "";
      char *fs3 = xlate_dna && stype==Combine ? 
          " (or its reverse)" : 
        stype==Combine ? 
          " (or its reverse complement)" : 
          "";
      char *fs4 = xlate_dna ? "," : ", and";		/* frame string 4 */
      char *fs5 = xlate_dna && stype==Combine ? 
	  ", and\n\t   o the protein translation of the match (or its reverse).\n" :
          xlate_dna ? 
	  ", and\n\t   o the protein translation of the match.\n" :
          ".\n";
      char *ptype = use_seq_p ? "SEQUENCE" : "POSITION";
      printf(
"\t- The positions and p-values of the non-overlapping motif occurrences\n"
"\t  are shown above the actual sequence for each of the high-scoring\n"
"\t  sequences from Section I.\n"
"\t- A motif occurrence is defined as a position in the sequence whose\n"
"\t  match to the motif has %s p-value less than %g as \n"
"\t  defined in Section II.\n"
"\t- For each sequence, the first line specifies the name of the sequence.\n"
"\t- The second (and possibly more) lines give a description of the \n"
"\t  sequence.\n"
"\t- Following the description line(s) is a line giving the length, \n"
"\t  combined p-value, and E-value of the sequence as defined in Section I.\n"
"\t- The next line reproduces the motif diagram from Section II.\n"
"\t- The entire sequence is printed on the following lines.\n"
"\t- Motif occurrences are indicated directly above their positions in the\n"
"\t  sequence on lines showing\n"
"\t   o the motif number%s of the occurrence%s,\n"
"\t   o the position p-value of the occurrence,\n"
"\t   o the best possible match to the motif%s%s\n"
"\t   o columns whose match to the motif has a positive score (indicated \n"
"\t     by a plus sign)%s",
      ptype, w_thresh, fs1, fs2, fs3, fs4, fs5);
      printf("%s\n", stars);
    } /* doc */
    copy_file(note_file, stdout);
    printf("\n%s\n\n", stars);
  } /* annotation table */

  /* display elapsed time */
  fflush(stdout);
  system("echo ''; echo CPU: `hostname`;");
  printf("Time %f secs.\n\n", myclock()/1E6);
} /* print_mast_results */

/**********************************************************************/
/*
	make_mast_tables

	Make the three MAST tables in temporary files.

	Returns the number of sequences less than e_max.
*/
/**********************************************************************/
static int make_mast_tables(
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of different strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  BOOLEAN use_seq_comp,			/* compensate for seq composition */
  BOOLEAN best_motifs, 			/* diagrams have only best motif */
  int kmotifs,	 			/* number of known motifs */
  MOTIF motif[],			/* data returned by read_motifs */
  char *motif_file,			/* motif file name */
  FILE *fdata,				/* the database file */
  int nseqs,				/* number of sequences scored */
  SORT_SEQ *seqs,			/* sortable array of seq/score records*/
  int n_hits,				/* size of seqs array */
  FILE* hit_file,			/* table of hits */
  FILE* diag_file,			/* table of motif diagrams */
  FILE* note_file,			/* table of annotated sequences */
  int rank,				/* rank of first result to print */
  int smax,				/* maximum sequences to print*/
  double e_max,				/* maximum E-value to print */
  double m_thresh,			/* maximum motif p-value to print */
  double w_thresh,	 		/* maximum motif p-value-- weak hits */
  BOOLEAN use_seq_p,			/* use sequence not position p-values */
  BOOLEAN hit_list,                     /* create hit list instead of diagram */
  LO *los[],				/* array of pointers to lo matrices */
  int nmotifs 				/* number motifs read */ 
)
{
  int i;
  int seqno; 				/* seq. number in sorted sequences */
  long length;				/* length of sample */
  char *sample_name;			/* name of sample */
  char *sequence;			/* sequence of sample */
  char *id;				/* identifier text for sample */
  int pseqs = 0;			/* number of sequences printed */

  for (seqno=0; seqno < n_hits; seqno++) {
    SORT_SEQ *seq = seqs + seqno;		/* next seq */
    TILING tiling;				/* tiling of sequence */
    int strand = seq->strand;			/* -1 neg, 0 both, +1 pos */
    double **pv = seq->pv;			/* pvalue tables */
    double pvalue = seq->Pvalue;		/* combined p-value of seq  */
    double evalue = nseqs * pvalue;		/* combined E-value of seq  */
    char name[1000];				/* name with gi|...| removed */

    /* done if E-value above threshold; rest have larger p-values */
    if (evalue > e_max) break;

    /* done if maximum number of sequences reached */
    if (smax && pseqs >= smax) break;

    /* skip if rank not reached */
    if (seqno+1 < rank) continue;
    pseqs++;

    fseek(fdata, seq->fp, 0);			/* go to start of sequence */
    read_sequence(fdata, &sample_name, &id, &sequence, &length);

#ifdef obsolete
    /* 
      Create a copy of the sequence name that has the "gi|...|" removed
      if there is more to the name.  Full name is used only in annotation
      section.
    */
    /* find amount to remove */
    if (strncmp(sample_name, "gi|", 3) == 0) { 	/* starts with "gi|" */
      char *s = strchr(sample_name+3, '|');	/* next "|" character */
      if (s != NULL) shift = (int) (s-sample_name+1);	/* shift amount */
    }
    /* copy the (truncated) name into the new location */
    for (i=shift; sample_name[i] && i-shift+1<1000; i++) {
      name[i-shift] = sample_name[i];
    }
    /* remove trailing "|" */
    if (i-shift+1<1000 && name[i-shift-1] == '|') i--;	
    name[i-shift] = '\0';
#endif
    strncpy(name, sample_name, 999);

    /*
      convert to reverse complement if negative strand of DNA
    */
    if (strand == -1) invcomp_dna(sequence, length);

    /* 
      score, tile and diagram the sequence with each of the motifs
    */
    tiling = score_tile_diagram(sequence, length, los, nmotifs, dna, stype, 
      xlate_dna, best_motifs, FALSE, pv, m_thresh, w_thresh, use_seq_p, hit_list);

    /* 
      Print the hit in the hits section. 
    */
    if (hit_file != NULL) {
      char *kp = !kmotifs ? "" : (seq->known_pos ? "+ " : "- ");  /* +/- */
      int nl = !kmotifs ? MMSN : MMSN - 2;		/* length of name */
      char *frame = xlate_dna ? best_frame(nmotifs, length, tiling) : "";
      char *st = stype==Separate ? (strand==-1 ? " -" : " +") : ""; 
      int il = (*st=='\0') ? MAXID : MAXID-2;		/* length of id */
      char *elipsis = strlen(id)>il ? "... " : "   ";
      if (*frame!='\0') il -= 2;
      fprintf(hit_file, "%s%-*.*s %-*.*s%3s%s%s  %8.2g %6ld\n", kp, nl, nl, 
        name, il, il, id, elipsis, st, frame, evalue, length);
    } /* print hit */

    /* 
      Print the motif diagram in the diagrams section.
    */
    if (diag_file != NULL){
      char hdr[80];						/* header */
      char *st = stype!=Separate ? "" : (strand==-1 ? " -" : " +"); 
      sprintf(hdr, "%-*.*s%s %8.2g  ", MMSN, MMSN, name, st, evalue);
      print_diagram(tiling.diagram, hdr, diag_file);
    }

    /* 
      Print annotated sequence in the annotation section.
    */
    if (note_file) print_annotation(dna, stype, xlate_dna, m_thresh, note_file,
      los, sequence, length, sample_name, id, strand, pvalue, evalue, tiling);

    /* 
      free space 
    */
    myfree(sample_name);
    myfree(id);
    myfree(sequence);
    free_tiling(tiling);
    if (use_seq_comp) {
      for (i=0; i<nmotifs; i++) myfree(seq->pv[i]);
      myfree(seq->pv);
    }
  } /* sequence */

  return(pseqs);
} /* make_mast_tables */

/**********************************************************************/
/*
	get_scores

	Calculate the score and p-value for each sequence in the
	database and sort them by p-value in ascending order. 

	Returns sorted array of sequence-score records.
*/
/**********************************************************************/
static SORT_SEQ *get_scores(
  BOOLEAN dna,		/* database is DNA */
  STYPE stype,		/* handling of DNA strands */
  BOOLEAN xlate_dna,	/* database is DNA and motifs protein */
  double f[],		/* null letter probability distribution */
  FILE *fdata,		/* the database */
  FILE *fsave,		/* good sequences if database on stdin */
  LO *los[],		/* array of pointers to log-odds matrices */
  int nmotifs,		/* number of log-odds matrices in los */
  int range, 		/* set entries in matrices to [0..range] */
  double **pv,		/* p-value tables for each motif */
  BOOLEAN sonly,	/* calculate p-value of observed spacings only */
  BOOLEAN lump,		/* combine spacings into one p-value */
  BOOLEAN use_seq_comp,	/* adjust E-value for actual sequence composition */
  double e_max,		/* maximum sequence E-value to print */
  int kmotifs,		/* number of known motifs*/
  MOTIF motif[],	/* data returned by read_motifs */
  BOOLEAN status,	/* show progress */
  int min_seqs,		/* lower bound on nseqs */
  int *nseqs,		/* total sequences in database */
  double *residues,	/* total number of residues in seqs */
  int *n_hits		/* sequences less than e_max */
)
{
  int i; 
  int imotif;
  long fp = 0;				/* pointer in datafile */
  char *sample_name;			/* name of sample */
  char *sequence;			/* sequence of sample */
  char *id;				/* identifier text for sample */
  long length;				/* length of the sequence */
  double *best_score = NULL;		/* best score per motif */
  int *best_loc = NULL;			/* best location per motif */
  SORT_SEQ *seqs = NULL;		/* sortable array of seq/score records*/
  double pvalue;			/* expected number of sequences >= x */
  SCORE **scores;			/* scores for each motif vs seq. pos. */
  int strand = (stype==Separate || stype==Norc) ? 1 : 0;	/* current */
  int alen = los[0]->alen;		/* length of motif alphabet */
  BOOLEAN saved = FALSE;		/* sequence already saved */

  *n_hits = *nseqs = *residues = 0;	/* actual and saved sequences */
  while (1) {

    /* 
      read next sequence unless doing - strand; break on EOF or error 
    */
    if (strand != -1) {					/* not - strand */
      fp = fsave ? ftell(fsave) : ftell(fdata);		/* save seq position */
      if (!read_sequence(fdata, &sample_name, &id, &sequence, &length)) break;
      saved = FALSE;					/* new sequence */
    }

    /*
      create space for best scores and locations per motif unless they exist
    */
    if (!best_score) Resize(best_score, nmotifs, double);
    if (!best_loc) Resize(best_loc, nmotifs, int);

    /* 
      update size of database 
    */
    if (stype==Separate && strand==1) {	/* treat each DNA sequence as two */
      (*nseqs)+=2;			/* number of sequences in database */
      *residues += 2*length;		/* number of residues in database */
    } else if (stype!=Separate) {	/* treat each sequence as one */
      (*nseqs)++;			/* number of sequences in database */
      *residues += length;		/* number of residues in database */
    }

    /*
      convert DNA to negative strand if just did the positive strand
    */
    if (strand == -1) invcomp_dna(sequence, length);

    /* 
      compute the motif scores
    */
    scores = score_sequence(stype, xlate_dna, sequence, length, nmotifs, los);

    /*
      find best scoring position for each motif for the sequence
    */
    for (imotif=0; imotif<nmotifs; imotif++) {	/* motif */
      long seq_pos;
      int ws = los[imotif]->ws;			/* width of motif in sequence */
      best_score[imotif] = LITTLE;
      best_loc[imotif] = 0;
      if (ws > length) continue;		/* skip if motif too long */
      for (seq_pos=0; seq_pos<=length-ws; seq_pos++) {	/* position in sequence */
	double s = scores[imotif][seq_pos].score;
	if (s > best_score[imotif]) {
	  best_score[imotif] = s;
	  best_loc[imotif] = seq_pos;
	}
      } /* position */
    } /* motif */
    free_2array(scores, nmotifs);		/* free scores array */

    /* 
      calculate the combined p-value for this sequence
    */
    pvalue = calc_p_value(sample_name, length, stype, nmotifs, los, best_score, 
      best_loc, range, pv, sonly, lump, norder, debug);

    /* 
      save the sequence if its E-value may be under threshold 
    */
    if (pvalue*(MAX(min_seqs, *nseqs)) < e_max) { 
      SORT_SEQ *seq;

      /* 
        create sequence record for sorting 
      */
      if (*n_hits % RCHUNK == 0) Resize(seqs, *n_hits+RCHUNK, SORT_SEQ);
      seq = &seqs[(*n_hits)++];		/* record for sorting */
      seq->id = sample_name;		/* name of sequence */
      seq->fp = fp;			/* start of sequence record */
      seq->length = length;		/* length of sequence */
      seq->Pvalue = pvalue;		/* p-value of product of p-values */
      seq->strand = strand;		/* strand of DNA or 0 if both */
      seq->score = best_score;		/* best score per motif */
      seq->loc = best_loc;		/* locations of best scores */
      seq->comp = use_seq_comp ? get_seq_comp(xlate_dna, sequence, alen) : NULL;
      seq->pv = pv;			/* pvalue tables for each motif */
      best_score = NULL;
      best_loc = NULL;			/* don't reuse the space! */

      /*
        write the sequence to the fsave file if reading standard input
      */
      if (fdata == stdin && !saved) {
        fprintf(fsave, ">%s %s\n%s\n", sample_name, id, sequence);
        saved = TRUE;			/* don't save this seq again */
      }

      /* 
	set flag if sequence is known positive and doing ROC 
      */
      if (kmotifs) {
	if (hash_lookup(sample_name, 1, motif[0].ht)) {
	  seq->known_pos = TRUE;
	} else if (hash_lookup(sample_name, 0, motif[0].ht)) {
	 /* skip sequences labeled col=0 in .tag file; "?/P" in Prosite */
	  (*n_hits)--;
	} else {
	  seq->known_pos = FALSE;
	}
      } /* ROC */

    } /* hit */

    /* 
      print progress report
    */
    if (status && *nseqs % 100 == 0) 
      fprintf(stderr, "\rsequences: %6d ", *nseqs);

    /*
      flip strand flag
    */
    if (stype==Separate) strand *= -1;

    /* 
     free up space 
    */
    if (strand != -1) {				/* don't free if first strand */
      if (!kmotifs) myfree(sample_name);	/* don't free if doing ROC */
      myfree(sequence);
      myfree(id);	
    }

  } /* read_sequence */

  /*
    recalculate the p-values based on actual sequence composition
  */
  for (i=0; use_seq_comp && i < (*n_hits); i++) {
    int j;
    double **pv = NULL;
    long length = seqs[i].length;
    SORT_SEQ *seq = seqs + i; 
    double *best_score = seq->score;
    int *best_loc = seq->loc;

    /* skip if E-value too big */
    if (seq->Pvalue * (*nseqs) > e_max) continue;	

    /*
      calculate p-values of all integer score values in range [0...w*RANGE] 
    */
    Resize(pv, nmotifs, double *);
    for (j=0; j<nmotifs; j++) {
      /* pv[j] = calc_cdf(los[j], range, seq->comp);*/
      pv[j] = calc_pssm_cdf(los[j]->w, los[j]->alen, range, los[j]->logodds, seq->comp);
    }

    /* 
      recalculate the combined p-value for this sequence using the actual
      sequence composition as the background model
    */
    seq->Pvalue = calc_p_value(NULL, length, stype, nmotifs, los, best_score, 
      best_loc, range, pv, sonly, lump, norder, debug);

#ifdef DEBUG
    if (isnan(seq->Pvalue)) {
      fprintf(stderr, "nan: i %d id %s pvalue %f\n", i, seq->id, seq->Pvalue); 
      for (j=0; j<alen; j++) fprintf(stderr, "%c %f\n", unhash(j), seq->comp[j]);
      abort();
    }
#endif

    /* 
      print progress report
    */
    if (status && i && i % 100 == 0) 
      fprintf(stderr, "\rrecalc p-value sequences: %6d ", i);

    /* 
      free pv space if E-value too big, otherwise
      save composition-based pv distribution
    */
    if (seq->Pvalue * (*nseqs) > e_max) {
      for (j=0; j<nmotifs; j++) myfree(pv[j]);
      myfree(pv);
    } else {
      seq->pv = pv;
    }
  } /* recalculate p-values */

  /* 
    sort the sequences by p-value ascending order 
  */
  qsort((char *)seqs, *n_hits, (int)sizeof(SORT_SEQ), p_compare);

  /* 
    bail if no sequences were read succesfully 
  */
  if (*nseqs == 0) {
    fprintf(stderr, "Quitting due to errors or empty database.\n"); 
    exit(1);
  }

  return seqs;
} /* get_scores */

/**********************************************************************/
/*
	get_motifs

	Read in the log-odds matrices (motifs) from standard input, 
        remove unused motifs and interpret the ordering and spacing diagram.
	Exits if there is an error in the input.

	Returns number of motifs in the (compacted) los array.
	Updates los as well as globals order and space.
	Updates s2b.
*/
/**********************************************************************/
static int get_motifs(
  BOOLEAN xlate_dna,		/* database is DNA and motifs protein */
  char *nlogodds,		/* name of log-odds matrix */
  char *alphabet,		/* alphabet of logodds matrices */
  char *blast_alphabet,		/* corresponding BLAST alphabet */
  int *p[MAXASCII],		/* alphabet permutation/substitution matrix */
  BOOLEAN shuffle,		/* shuffle the motif columns */
  double min_ev,		/* minimum E-value of motifs to use */
  int umotifs, 			/* number motifs to be used */
  int mlist[MAXLO],		/* list of motifs given in -m option */
  BOOLEAN motifs[MAXLO],	/* list of motifs to output */
  LO *los[MAXLO],		/* array of logodds matrices */
  double *f,			/* array of null letter probabilities */
  int range 			/* set logodds matrices to [0..range] */
)
{
  int i, imotif, cnt;
  int nmotifs;			/* number of log-odds matrices in los */

  /* 
    read in log-odds matrices 
  */
  nmotifs = read_log_odds(xlate_dna, nlogodds, alphabet, blast_alphabet,
    p, range, los, f);

  /*
    check that at least one motif was read
  */
  if (nmotifs == 0) {
    fprintf(stderr, "No scoring matrices found.\n");
    exit(1); 
  }

  /* 
    clear the list of flags showing motifs to use
  */
  for (i=0; i<MAXLO; i++) motifs[i] = FALSE;

  /* 
   set flags of motifs to use; check for valid motif numbers 
  */
  if (umotifs > 0) {			/* using specified motifs */
    for (i=0; i<umotifs; i++) {
      int m = mlist[i];
      if (m < 1 || m > nmotifs) {
        fprintf(stderr, "Motif %d in -m option not in legal range: 1 to %d.\n",
          m, nmotifs); 
        exit(1);
      }
      motifs[m-1] = TRUE;
    }
  } else {				/* using motifs with E-value < min_ev */
    for (i=0; i<nmotifs; i++) {
      motifs[i] = (los[i]->ev < min_ev);
      umotifs++;
    }
  }

  /*
    flag motifs with no information
  */
  for (i=0; i<nmotifs; i++) {
    if (motifs[i] && !los[i]->scale) {		/* motif has no info */
      motifs[i] = FALSE;
      umotifs--;
    }
  }

  /*
    check that valid motifs remain
  */
  if (!umotifs) {
    fprintf(stderr, "No scoring matrices contained any information.\n");
    exit(1); 
  }

  /* 
    Parse the ordering and spacing diagram.
    Variables norder, order and space are globals defined in diagram.h 
  */

  if (diagram == NULL) {
   norder = 0;					/* no diagram */
  } else {					/* have a diagram */
    BOOLEAN fail = FALSE;			/* ok */
    BOOLEAN used[MAXLO];			/* for checking use once */

    dptr = 0;					/* diag read pointer */
    norder = 0;					/* lnumber of motifs in diag */
    for (i=0; i<MAXLO; i++) space[i] = 0;	/* set all spacers to zero */ 
    if (yyparse()) exit(1);			/* parse diagram; sets globals
						   norder, order, space */
    /* check that no unused motifs are in motif diagram */
    for (i=0; i<norder; i++) {
      int m = order[i];				/* motif number (plus 1) */
      if (m < 1 				/* number too small */
          || m > nmotifs			/* number too large */
          || !motifs[m-1]			/* motif not being used */
      ) {
	fprintf(stderr, "Unknown or unused motif %d given in motif diagram.\n",
          m);
        fail = TRUE;
      }
    }

    /* check that each motif only used once in diagram */
    for (i=0; i<MAXLO; i++) used[i] = FALSE; 
    for (i=0; i<norder; i++) {
      int m = order[i];
      if (used[m]) {
	fprintf(stderr, 
          "Motif %d used more than once in motif diagram:\n  %s\n", m, diagram);
	fail = TRUE;
      } else {
	used[m] = TRUE;
      }
    }
    
    /* quit if there was a problem in the diagram */
    if (fail) exit(1);

    /* change format for spacers to be relative to start of prior motif */
    space[0] = -1;				/* ignore first spacer */
    for (i=1; i<norder; i++) {
      int m = order[i-1];			/* previous motif in diagram  */
      space[i] += los[m-1]->w; 			/* internal format is -1 */
    }
  }
  
  /* 
    remove any motifs that are not to be used 
  */
  for (imotif=0, cnt=0; imotif<nmotifs; imotif++) {	/* motif */
    if (motifs[imotif]) {
      los[cnt++] = los[imotif];			/* logodds matrix */
    }
  } /* motif */
  nmotifs = cnt;				/* number of used motifs */

  /* 
    convert spacing diagram to internal motif names 
  */
  for (i=0; i<norder; i++) {
    for (imotif=0; imotif < nmotifs; imotif++) {
      if (order[i] == los[imotif]->name) break;	/* found motif name */
    }
    order[i] = imotif;				/* internal name for motif */
  }

  /* shuffle the columns of each motif matrix if requested */
  if (shuffle) shuffle_cols(los, nmotifs);

  /* compute the pairwise motif correlations (similarities) */
  motif_corr(nmotifs, los);
 
  /* print motif correlations */
  if (debug) {
    for (i=1; i<nmotifs; i++) {
      int j;
      printf("SCAL %2d", i+1);
      for (j=0; j<i; j++) {
        printf(" %8.2f", los[i]->corr[j]);
      }
      printf("\n");
    }
  }

  /* return number of motifs remaining in (compacted) los array */
  return nmotifs;
} /* get_motifs */

/**********************************************************************/
/*
	print_annotation

	Print the annotation section of MAST output.

	Format:
		<Sequence Description line>+
		<length> <combine p-value> <E-value>
		<block diagram>
		[[<strand>]<motif number>[<frame>]
		 <position p-value>
		 <best possible match>
		 <per-letter positive score markers>
		 [protein translation of motif if xlate_dna]
		 <sequence>
		]+

        Hits:
		strong 		p-value < thresh 
		weak		otherwise
*/
/**********************************************************************/
static void print_annotation(
  BOOLEAN dna,				/* database is DNA */
  STYPE stype,				/* handling of DNA strands */
  BOOLEAN xlate_dna,			/* database is DNA and motifs protein */
  double thresh,			/* threshold for strong hits */
  FILE *note_file,			/* table of annotated sequences */
  LO *los[],				/* array of pointers to lo matrices */
  char *sequence,			/* sequence of sample */
  long length,				/* length of sample */
  char *sample_name,			/* name of sample */
  char *id,				/* identifier text for sample */
  int strand,				/* strand of DNA or 0 if protein */
  double pvalue,			/* combined p-value of sequence */
  double evalue,			/* combined E-value of sequence */
  TILING tiling				/* tiling and diagram of sequence */
)
{
  long i, j, k;
  int gb;				/* no good break yet */
  int nw;				/* width of sequence number */
  int lw;				/* amount of seq per line */
  char *s;				/* pointer to start of id */
  int idlines = 0;			/* number of id lines printed */
  char *protalph = PROTEINB;		/* protein alphabet */
  char *bstring = NULL;			/* hit motif number string */
  char *pstring = NULL;			/* hit p-value string */
  char *cstring = NULL;			/* hit consensus string */
  char *mstring = NULL;			/* hit match string */
  char *xstring = NULL;			/* hit match string */
  char *sstr = stype!=Separate ? "" : 
    (strand==-1 ? " (- strand)" : " (+ strand)");
  int inc = (xlate_dna ? 3 : 1);	/* increment in sequence */

  /*
    Allocate space for annotation strings.
  */
  Resize(bstring, 2*length, char);
  Resize(pstring, 2*length, char);
  Resize(cstring, 2*length, char);
  Resize(mstring, 2*length, char);
  Resize(xstring, 2*length, char);

  /* 
    print the description of the sequence
  */
  fprintf(note_file, "\n\n%s%s\n", sample_name, sstr);
  for (s=id, gb=0; ; s+=(gb+1), gb=0) {		/* more to print */
    /* find the next good breaking point */
    for (i=0; i<PAGEWIDTH-2; i++) {
      if (s[i]==' ' || !s[i]) gb = i;		/* good break: blank or EOS */
      if (!s[i]) break; 			/* end of sequence reached */
    }
    if (gb == 0) gb = i-1;			/* no good break; back up 1 */
    fprintf(note_file, "  %-.*s\n", gb+1, s);	/* print to break inclusive */ 
    if (!s[i]) break;				/* done with sequence */
    if (++idlines == MAXIDLINES) {fprintf(note_file,"  ...\n"); break;}
  }
	
  /* 
    print the length, combined p-value and E-value 
  */
  fprintf(note_file, 
    "  LENGTH = %ld  COMBINED P-VALUE = %8.2e  E-VALUE = %8.2g\n",
    length, pvalue, evalue);

  /* 
    print the motif diagram 
  */
  print_diagram(tiling.diagram, "  DIAGRAM: ", note_file);
  putc('\n', note_file);			/* blank line */

  /* 
    print the sequence and its motif annotation 
  */

  /* create lines holding the motif numbers, p-values, consensus sequences
     and matches for each hit 
  */
  for (i=0; i < length; i++) {			/* position in sequence */
    int m = abs(tiling.hits[i]) - 1;		/* motif number of hit */

    if (m >= 0) {				/* motif at position i */
      LO *lo = los[m];
      int w = lo->w;				/* width of motif */
      int ws = lo->ws;				/* width of motif in sequence */
      LOGODDS logodds = lo->logodds;		/* logodds matrix */
      double scale = lo->scale;			/* scale for convert 2 bits */
      double offset = lo->offset;		/* offset for convert 2 bits */
      char *p;					/* positive match character */
      char *s = sequence+i;			/* start of hit */
      char mns[80];				/* motif number string */ 
      BOOLEAN ic = tiling.hits[i] < 0;		/* hit on - strand */
      char *cons = ic ? los[m]->best_icseq : los[m]->best_seq;	/* consensus */
      char *pfmt = ws>6 ? "%-*.1e" : "%-*.0e";	/* pvalue format */
      char *mfmt = ic ? "%*.*s" : "%-*.*s";	/* match format */
      char *cfmt = xlate_dna ? (ic ? "..%c" : "%c..") : "%c";	/* cons. fmt. */
      char *strand = stype==Combine ? (ic ? "-" : "+") : ""; 	/* hit strnd */
      int frame = xlate_dna ? i%3 + 1 : 0;	/* frame if translating DNA */

      /* make the block diagram */
      make_block(los[m]->name, strand, frame, thresh, tiling.pvalues[i], 
        FALSE, mns);

      /* align things in separate strings */
      sprintf(bstring+i, "%-*.*s", ws, ws, mns);	/* motif numbes */
      sprintf(pstring+i, pfmt, ws, tiling.pvalues[i]);	/* p-value */

      /* consensus, positive matches (, translation to protein) */
      for (j=0, k=i; j<w; j++, k+=inc, s+=inc) {
        int cx = chash(xlate_dna, ic, s);
        double score = ic ? logodds(w-j-1, cx) : logodds(j, cx); 
        sprintf(cstring+k, cfmt, cons[j]);		/* consensus */
	p = (scaled_to_bit(score, 1, scale, offset) > 0) ? "+" : " ";
        sprintf(mstring+k, mfmt, inc, inc, p);		/* positive match */
        if (xlate_dna) sprintf(xstring+k, cfmt, protalph[cx]);	/* xlated */
      }
      i += ws - 1;
    } else {				/* no motif at position i */
      bstring[i] = pstring[i] = cstring[i] = mstring[i] = xstring[i] = ' ';
    } /* motif at position i */
  } /* position in sequence */
  bstring[i] = pstring[i] = cstring[i] = mstring[i] = xstring[i] = '\0';

  /* 
     Print the annotation lines for this sequence. 
     If the sequence line contains a motif occurrence print 
	    1) blank line (if not first line)
	    2) motif number line
	    3) consensus line
	    4) match line
	    5) [DNAB to PROTEINB translation of motif]
	    6) sequence line
     otherwise print only
	    1) blank line (if last line had motif)
	    2) sequence line
  */
  nw =  MAX(4, 1+(int)(log(length)/log(10.0)));	/* width of position number */
  lw = PAGEWIDTH - (nw + 1);			/* amount of sequence to print*/
  for (j=0; j < length; j += lw) {		/* position in sequence */
    int p = (strand == -1) ? length-j : j+1;	/* number - strand backwards */
    for (k=0; k<lw; k++) if (cstring[j+k] != ' ') break; 
    if (k < lw && cstring[j+k] != '\0') {	/* line contains a motif */
      char *bs = bstring+j;			/* start of block line */	
      char *ps = pstring+j;			/* start of pvalues line */	
      char *cs = cstring+j;			/* start of consensus line */	
      char *ms = mstring+j;			/* start of matches line */	
      char *xs = xstring+j;			/* start of xlated motif line */
      int bw, pw, cw, mw, xw;			/* width of various lines */
      /* find lengths of lines */
      bw = pw = cw = mw = xw = MIN(lw, length-j);
      for (; bw>0 && bs[bw-1]==' '; bw--);
      for (; pw>0 && ps[pw-1]==' '; pw--);
      for (; cw>0 && cs[cw-1]==' '; cw--);
      for (; mw>0 && ms[mw-1]==' '; mw--);
      if (xlate_dna) for (; xw>0 && xs[xw-1]==' '; xw--);
      /* print lines */
      if (j != 0) putc('\n', note_file);			/* blank */
      fprintf(note_file, "%*.*s %*.*s\n", nw,nw,"", bw,bw,bs);	/* motifs */
      fprintf(note_file, "%*.*s %*.*s\n", nw,nw,"", pw,pw,ps);	/* pvalues */
      fprintf(note_file, "%*.*s %*.*s\n", nw,nw,"", cw,cw,cs);	/* consensus */
      fprintf(note_file, "%*.*s %*.*s\n", nw,nw,"", mw,mw,ms);	/* matches */
      if (xlate_dna) 
        fprintf(note_file, "%*.*s %*.*s\n", nw,nw,"", xw,xw,xs);/* xlt'd mtf */
      fprintf(note_file, "%-*d %.*s\n", nw, p, lw, sequence+j);	/* sequence */
    } /* line contains motif */
  } /* position in sequence */

  /*
   Free space.
  */
  myfree(bstring);
  myfree(pstring);
  myfree(cstring);
  myfree(mstring);
  myfree(xstring);

} /* print_annotation */

/**********************************************************************/
/*
	best_frame

	Find the predominant frame of the matches to the motifs for
	a sequence.  The best frame is the frame with the minimum
	product of p-values.

	Returns a string "a ", "b " or "c ".
*/
/**********************************************************************/
static char *best_frame(
  int nmotifs, 				/* number of motifs */
  long length,				/* length of sequence */
  TILING tiling				/* tiling of sequence */
)
{
  int i, m;
  long j;
  double best_frame_p = 1;		/* p-value of best frame */
  int best = 0;				/* best frame */

  /* find the product of p-values for each frame */
  for (i=0; i<3; i++) {				/* frame */
    double best_p[MAXLO];			/* best p-values for motifs */
    double prod_p = 1;				/* product of p-values */

    /* clear array of best p-values for each motif */
    for (j=0; j<nmotifs; j++) best_p[j] = 1;	/* best p-value for motif j */

    /* find best p-value for each motif in this frame */
    for (j=i; j<length; j+=3) {			/* position in sequence */
      if ((m = abs(tiling.hits[j])-1) >= 0) 
        best_p[m] = MIN(best_p[m], tiling.pvalues[j]);
    }

    /* form the product of p-values for this frame */
    for (j=0; j<nmotifs; j++) prod_p *= best_p[j];

    /* update best frame */
    if (prod_p < best_frame_p) {
      best_frame_p = prod_p;
      best = i;
    }
  } /* frame */
 
  return (best == 0 ? "a " : (best == 1 ? "b " : "c "));
} /* best_frame */

/************************************************************************/
/*
	get_bad_motifs

	Get list of motifs to remove so that no pairs of motifs will
	have excessive correlation.

	Returns number of bad motifs and the list of bad motifs.
*/
/************************************************************************/
static int get_bad_motifs(
  double maxcorr,			/* maximum allowed corellation */
  int nmotifs,				/* number of motifs */ 
  LO *los[],				/* array of logodds matrices */
  int bad_list[] 			/* list of highly correlated motifs */
)
{
  int i, j;
  int nbad = 0;				/* number of highly correlated motifs */

  for (i=1; i<nmotifs; i++) {		/* from motif */
    BOOLEAN ibad = FALSE;		/* motif highly correlated */
    for (j=0; j<i; j++) {		/* to motif */
      if (!ibad && los[i]->corr[j] > maxcorr) {
	ibad = TRUE;
	bad_list[nbad++] = los[i]->name;
      }
    } /* to motif */
  } /* from motif */

  return nbad;
} /* get_bad_motifs */

