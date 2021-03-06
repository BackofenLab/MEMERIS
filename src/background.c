/*
 * $Id: background.c,v 1.1.1.1 2005/07/25 23:49:10 nadya Exp $
 * 
 * $Log: background.c,v $
 * Revision 1.1.1.1  2005/07/25 23:49:10  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/****************************************************************************
*                                                                           *
*       MEME++                                                              *
*       Copyright 1994-1999, The Regents of the University of California    *
*       Author: Timothy L. Bailey                                           *
*                                                                           *
****************************************************************************/

#include "macros.h"
#include "background.h"
#include "hash_alph.h"

/* index from ASCII to integer and back */
static int a2i_index[MAXASCII];
static int i2a_index[MAXASCII];
static int index_alen;				/* length of alphabet */
#define a2i(a) a2i_index[(int)(a)]
#define i2a(i) i2a_index[(int)(i)]

static int setup_index(
  char *alpha					/* the alphabet */
);

static int s2i(
  char *string					/* the string to index */
);

static char *check_prob(
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  double  *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
); 

static void average_rc(
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  double *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
);

static void add_x_tuples(
  double *a_p,					/* tuple->prob assoc. array */
  char *ntuple,					/* non-X-containing tuple */
  int ntuplew, 					/* ntuple width */
  int ntuplei, 					/* ntuple index */
  char *xtuple,					/* X-containing tuple */
  int xtuplew, 					/* xtuple width */
  int xcnt					/* count of X's in xtuple */
);

static void print_prob(
  double *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
);
 
static double *get_cond_prob (
  double *a_p, 					/* tuple->prob array */
  int n						/* length of array */
);

/***************************************************************************/
/*
	read_markov_model

	Read a file of tuple probabilities defining a Markov model.  
	Check that all required tuples are present.  

	File must have format:
		[<tuple> <p>]+
	and may have comment lines beginning with "#" in column 1.

	If the file is NULL and freq is not, uses the 0-order
	frequencies in freq as the model, after adding X if requested.

	Sets the order of the model.

	Returns array: 
		cp[s2i(wa)] = Pr(a | w)
	where "a" is a character and "w" is a string.
	The 0-order probabilities are in positions 0..alength-1 of the
	array.

*/
/***************************************************************************/
extern double *read_markov_model( 
  char *pfile, 					/* name of probability file */
  double *freq,					/* letter frequencies */
  char *alpha,					/* alphabet expected */
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  BOOLEAN rc,					/* average reverse complements*/
  int *order					/* order of model read */
) 
{
  int i;					/* index into array */
  double a_p[MAX_BACK_SIZE];			/* tuple-prob array */
  double *a_cp=NULL; 				/* conditional prob. array */
  FILE *pfilep;					/* file pointer to file */
  char *line=NULL;				/* line buffer */
  char **fields=NULL;				/* fields of line */
  int nfields;					/* number of fields in line */
  int line_no=0;				/* line number */
  char *tuple;					/* the tuple */
  double p;					/* the probability */
  int maxw=0;					/* maximum tuple width */
  int alen=strlen(alpha);			/* length of alphabet */
  int ntuples;					/* number of tuples */

  /* check input */
  if (!pfile && !freq) {
    fprintf(stderr, "read_markov_model error: specify pfile or freq\n");
    exit(1);
  }

  /* add 'X 'to the alphabet if requested */
  if (add_x) {
    char *tmp = NULL;
    Resize(tmp, alen+2, char);
    strcpy(tmp, alpha);
    tmp[alen] = 'X'; tmp[alen+1] = '\0';
    alpha = tmp;
    alen++; 
  }

  /* setup the mapping from ascii to integer and back */
  setup_index(alpha);

  /* use the frequencies if given */
  if (freq) {					/* frequencies given */
    Resize(a_cp, alen, double);
    for (i=0; i<alen-add_x; i++) RND(freq[i], 8, a_cp[i]);
    if (add_x) a_cp[i] = 1.0;			/* Pr(X) */
    /* averate reverse complement probabilities together if requested */
    if (rc) average_rc(add_x, a_cp, 1, "", 0, alpha); 
    return(a_cp);
  }

  /* initialize probability array */
  for (i=0; i<MAX_BACK_SIZE; i++) a_p[i] = -1;

  /* read in the probabilities and save indexed by uppercase tuple name */
  if (!(pfilep = fopen(pfile, "r"))) {
    fprintf(stderr, "Unable to open file %s for reading.\n", pfile);
    exit(1);
  }

  /*fprintf(stderr, "Reading background probabilities...\n");*/
  while (1) {					/* read file */
    int len, index;
    line_no++;
    Getline(pfilep, line, len);			/* read next line */
    if (!line) break;				/* at EOF */
    if (line[0] == '#') continue;		/* skip comments */
    Split(line, fields, nfields);		/* get tuple and prob */
    if (nfields != 2) {
      fprintf(stderr, 
        "Formatting error in file %s line %d: %s\n", pfile, line_no, line);
      exit(1);
    }
    tuple = fields[0];
    p = atof(fields[1]);
    if (p<0 || p>1) {
      fprintf(stderr, "Illegal probability in file %s line %d: %s\n", 
        pfile, line_no, line);
    }
    len = strlen(tuple);
    maxw = MAX(len, maxw);
    index = s2i(tuple);
    if (index < 0) {
      fprintf(stderr, "Illegal character in word `%s' in file %s line %d: %s\n",
        tuple, pfile, line_no, line);
      exit(1);
    }
    if (index >= MAX_BACK_SIZE) {
      for (i=1, ntuples=0; i<=maxw; i++) ntuples+= pow(alen, i);
      fprintf(stderr, "Background model too large.  Use smaller model or increase \nMAX_BACK_SIZE to at least %d in background.h and recompile.\n", ntuples);
      exit(1);
    }
    a_p[index] = p;				/* store probability */
  }
  fclose(pfilep);

  /* check that all necessary probabilities are defined */
  /*fprintf(stderr, 
    "Checking that all necessary probabilities are defined...\n");
  */
  tuple = check_prob(add_x, a_p, maxw, "", 0, alpha); 
  if (tuple) { 
    fprintf(stderr, "File %s gives no probability for %s.\n", pfile, 
      tuple);
    exit(1);
  }

  /* averate reverse complement probabilities together if requested */
  if (rc) average_rc(add_x, a_p, maxw, "", 0, alpha); 

  *order = maxw - 1;				/* order of Markov model */

  /* get conditional probabilities */
  for (i=1, ntuples=0; i<=maxw; i++) ntuples+= pow(alen, i);
  a_cp = get_cond_prob(a_p, ntuples);

  /* print the probabilities */
#ifdef DEBUG
  print_prob(a_cp, maxw, "", 0, alpha);
#endif

  return(a_cp);					/* return conditionals */
} /* read_markov_model */

/***************************************************************************/
/* 
 	check_prob
 
 	Recursively check that all probabilities are defined for all tuples
	up to width n.  Add their probabilities to all the X-containing
	tuples they match.
 
	Returns NULL if OK, otherwise the missing tuple. 
*/
/***************************************************************************/
static char *check_prob(
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  double *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
) 
{
  int i;
  char *t = NULL, *x = NULL; 			/* tuple and x-tuple */
  char *missing;				/* missing tuple */
  int ti;					/* index of tuple */

  if (n==0) return(NULL);			/* everything is OK */

  Resize(t, tuplew+2, char);
  Resize(x, tuplew+2, char);
  
  for(i=0; alpha[i+add_x]; i++) {		/* ignore last letter (X) */
    /* append letter to tuple */
    strcpy(t, tuple); t[tuplew] = alpha[i]; t[tuplew+1] = '\0';
    /* check that tuple exists */
    ti = s2i(t);				/* index of tuple */
    if (a_p[ti] < 0) return(t);
    /* add the current tuple's probability to matching X-containing tuples */
    if (add_x) add_x_tuples(a_p, t, tuplew+1, ti, x, 0, 0);
    /* recur */
    missing = check_prob(add_x, a_p, n-1, t, tuplew+1, alpha);
    if (missing) return(missing);
  } /* letter */
  myfree(t);
  if (add_x) myfree(x);

  return(NULL);					/* all found */
} /* check_prob */

/***************************************************************************/
/* 
 	average_rc
 
 	Recursively average the probabilities of reverse complement tuples.
 
*/
/***************************************************************************/
static void average_rc(
  BOOLEAN add_x,				/* add x-tuples if TRUE */
  double *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
)
{
  int i, j;
  char *t = NULL, *rc = NULL; 			/* tuple and rc-tuple */
  int ti, rci;					/* index of tuple */

  if (n==0) return;				/* everything is OK */

  Resize(t, tuplew+2, char);
  Resize(rc, tuplew+2, char);
  
  for(i=0; alpha[i+add_x]; i++) {		/* ignore last letter (X) */
    /* append letter to tuple */
    strcpy(t, tuple); t[tuplew] = alpha[i]; t[tuplew+1] = '\0';
    /* make reverse complement of tuple */
    for (j=0; j<=tuplew; j++) rc[j] = comp_dna(t[tuplew-j]);
    rc[tuplew+1] = '\0';
    /* get tuple and rc indices */
    ti = s2i(t);				/* index of tuple */
    rci = s2i(rc);				/* index of reverse complement*/
    /* average their probabilites */
    a_p[ti] = a_p[rci] = (a_p[ti] + a_p[rci]) / 2.0;
    /* recur */
    average_rc(add_x, a_p, n-1, t, tuplew+1, alpha);
  } /* letter */

  myfree(t);
  myfree(rc);

  return;
} /* average_rc */


/***************************************************************************/
/* 
 	print_prob
 
 	Recursively print the probabilities for all tuples
	up to width n.  
*/
/***************************************************************************/
static void print_prob(
  double *a_p,					/* tuple->prob assoc. array */
  int n,					/* widest tuple */
  char *tuple,					/* call with "" */
  int tuplew,					/* tuple width; call with 0 */
  char *alpha 					/* alphabet */
) 
{
  int i;
  char *t = NULL; 				/* tuple */

  if (n==0) return;				/* everything is OK */

  Resize(t, tuplew+2, char);
  
  for(i=0; alpha[i]; i++) {			/* include last letter (X) */
    /* append letter to tuple */
    strcpy(t, tuple); t[tuplew] = alpha[i]; t[tuplew+1] = '\0';
    /* print tuple and probability */
    printf("%10s %8.4f\n", t, a_p[s2i(t)]);
    /* recur */
    print_prob(a_p, n-1, t, tuplew+1, alpha);
  } /* letter */
  myfree(t);

} /* print_prob */

/***************************************************************************/
/*
	add_x_tuples

	Compute the probability of all X-containing tuples by
	adding the probability of the current non-X-containing tuple
	to each X-containing tuple it matches.  When this has
	been done for all non-X-containing tuples, the X-containing
	tuple probabilities will be correct.
*/
/***************************************************************************/
static void add_x_tuples(
  double *a_p,					/* tuple->prob assoc. array */
  char *ntuple,					/* non-X-containing tuple */
  int ntuplew, 					/* ntuple width */
  int ntuplei, 					/* ntuple index */
  char *xtuple,					/* X-containing tuple */
  int xtuplew, 					/* xtuple width */
  int xcnt					/* count of X's in xtuple */
)
{
  if (xtuplew == ntuplew) {			/* add prob. of non-X-tuple */
    int xtuplei;
    xtuple[xtuplew] = '\0';			/* end tuple */
    xtuplei = s2i(xtuple);
    if (a_p[xtuplei] == -1) a_p[xtuplei] = 0;	/* first time for x-tuple */
    a_p[xtuplei] += a_p[ntuplei];
  } else {					/* continue building x-tuple */
    xtuple[xtuplew] = 'X';			/* add X to xtuple */
    add_x_tuples(a_p, ntuple, ntuplew, ntuplei, xtuple, xtuplew+1, xcnt+1);
    if (xcnt || xtuplew < ntuplew-1) {		/* OK to add non-X to xtuple */
      xtuple[xtuplew] = ntuple[xtuplew];	/* add X to xtuple */
      add_x_tuples(a_p, ntuple, ntuplew, ntuplei, xtuple, xtuplew+1, xcnt);
    }
  }
} /* add_x_tuples */

/***************************************************************************/
/*
	get_cond_prob

	Compute 
		Pr(a | w)
	where "w" is a word and "a" is a character and store in an
	array with index "s2i(wa)".
*/
/***************************************************************************/
static double *get_cond_prob (
  double *a_p, 					/* tuple->prob array */
  int n						/* length of array */
)
{
  int wa, w;
  double *a_cp=NULL;

  Resize(a_cp, n, double);			/* create array */
  
  for (wa=0; wa<n; wa++) {			/* tuple "wa" */
    if (wa < index_alen) {			/* w = "" */
      a_cp[wa] = a_p[wa];			/* Pr(a) */
    } else {
      w = (wa/index_alen)-1;			/* get prefix w */
      a_cp[wa] = MIN(a_p[wa]/a_p[w], 1.0);	/* Pr(a | w) = Pr(wa)/Pr(w) */
    }
  }

  return(a_cp);
} /* get_cond_prob */

/***************************************************************************/
/*
 	log_cum_back

	Get log of the probability of each position in a string given a 
	Markov background model.

	Returns the cumulative background as an array:
		logcumback_i = 0, i=0
		             = log Pr(s_{0,i-1} | H_0), otherwise.
	and the total (log) cumulative background probability.

	The probability of any length-w substring starting at position i in the
	string (in the context of the string) can then be computed as:
		last_p = i+w-1;
		log_p = logcumback[last_p+1] - logcumback[i];

	The background model, H_0, is a Markov model defined by the
	order, n, and the conditional probabilities, a_cp, where
		a_cp[s2i(wa)] = Pr(a | w).

*/
/***************************************************************************/
extern double log_cum_back(
  char *seq,					/* the sequence */
  double *a_cp,					/* the conditional probs */
  int order,					/* the order of the model */
  double *logcumback				/* the result */ 
)
{
  int i;
  int len = strlen(seq);			/* length of sequence */

  /* compute probability for each position in string */
  for (i=0, logcumback[0]=0; i<len; i++) {	/* position of "a" */
    int lp = MIN(i+1, order+1);			/* substring length */
    char *wa;
    char *wa_end = &seq[i+1];
    char save;
    double log_pwa;				/* log Pr(a | w) */
#ifdef dumb
    Substr(wa, seq, i-lp+1, lp);		/* get wa */
#endif
    wa = &seq[i-lp+1];				/* next wa */
    save = *wa_end;				/* save character after wa */
    *wa_end = '\0';				/* make wa a string */
    log_pwa = log(a_cp[s2i(wa)]);		/* log Pr(a | w) */
    logcumback[i+1] = logcumback[i] + log_pwa;
    *wa_end = save;				/* restore seq */
/*
    myfree(wa);
*/
  }
  return(logcumback[i]);			/* total cumulative prob. */
} /* log_cum_back */

/***************************************************************************/
/*
	setup_index

	Setup the index from alphabet to integer and back.
	Upper and lower case index to same position; position indexes to
	upper case only.

	Returns length of alphabet.
*/
/***************************************************************************/
static int setup_index(
  char *alpha					/* the alphabet */
)
{
  int i, a;

  /* flag unused letters */
  for (i=0; i<MAXASCII; i++) a2i(i) = -1;	/* illegal letters */

  /* set up the hashing and unhashing indices */
  for (i = 0; (a = alpha[i]); i++) {
    a = islower(a) ? toupper(a) : a;		/* convert to uppercase */
    a2i(a) = i;					/* position in alphabet */
    i2a(i) = a;
    if (isupper(a)) a2i(tolower(a)) = i;	/* set lowercase as well */
  }
  index_alen = i;				/* set global */

  return(i);
} /* setup_index */

/***************************************************************************/
/*
	s2i

	Convert a string to an array index.

	Returns the index or -1 if the string contains illegal characters.
*/
/***************************************************************************/
static int s2i(
  char *string					/* the string to index */
)
{
  int i;
  char a;
  int index=0;					/* array index */

  /* complete the computation */
  while ((a = *(string++))) {
    i = a2i(a);					/* first charcter */
    if (i < 0) return(-1);
    index = index_alen*index + (i+1);
  }

  return(index-1);
} /* s2i */

/*
	To compile:
	gcc background.c -o background -I INCLUDE -lm -DSO
*/
#ifdef SO
main(
  int argc,
  char **argv
) 
{
  int i;
  char *pfile = argv[1];
  double logcumback[1000];
  char *string = "TTTTTTTTTTAAAAAAAAA";
  int order;
  int len = strlen(string);
  double *a_cp = read_markov_model(pfile, "ACGT", TRUE, &order);

  printf("order = %d\n", order);
  log_cum_back(string, a_cp, order, logcumback);
  for (i=0; i<len; i++) printf("%c %f %f\n", string[i], exp(logcumback[i]), 
    exp(Log_back(logcumback,i,6)) ); 
  printf("All done\n");

  return(0);
} /* main */
#endif /* SO */
