/*
 * $Id: read_seq_file.c,v 1.2.4.2 2006/01/31 20:30:46 nadya Exp $
 * 
 * $Log: read_seq_file.c,v $
 * Revision 1.2.4.2  2006/01/31 20:30:46  nadya
 * Attemp to fix control-m in description
 *
 * Revision 1.2.4.1  2006/01/31 19:11:54  nadya
 * add '\r' character to look for when delimiting strings
 * this is a fix for cygwin newline representation
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:24:50  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/***********************************************************************
*								       *
*	MEME++							       *
*	Copyright 1994, The Regents of the University of California    *
*	Author: Timothy L. Bailey				       *
*								       *
***********************************************************************/

/* 9-30-99 tlb; add sequence to SAMPLE */
/* 6-29-99 tlb; add resic */
/*
	Supports FASTA and SWISSPROT formats.

Note:
	Assumes SWISSPROT format if "ID" appears at beginning of first line.
	Ignores duplicate (same <id>) sequences.

All formats:
	<id>		:= <alpha>+
			   The *unique* name of the sequence; may
			   be no longer than 2*MSN characters; may *not*
			   contain whitespace!

			   NOTE: if <id> is "WEIGHTS", the <description>
			   field is parsed as a string of sequence weights 
			   and any <sequence> field (if any) is ignored.
			   All weights are assigned in order to the 
			   sequences in the file. If there are more sequences
			   than weights, the remainder are given weight one.
			   Weights must be greater than zero and less than
			   or equal to one.  Weights may be specified by
			   more than one "WEIGHT" entry which may appear
			   anywhere in the file.

	<description>	:= <alpha>+
			   A description of the sequence;
		           may contain whitespace.

	<sequence>	:= <alpha>+
			   The DNA or protein sequence data; 
                           may contain whitespace.

	<text>		:= <alpha>*
			   Text which is ignored.


Pearson FASTA format:
	INPUT DATA FILE FORMAT:
		[
		<header>
		<sequence>+
		]*

	<header>	:= ">"<id>" "<description>
			    Everything up to the first space is 
			    assumed to be the unique name of the sequence.
	
SWISSPROT Format
	[
	<ID>
	<DE>+
	<SQ>
	<sequence>+
	//
	]*

	<ID>		:= "ID   "<id>" "<text>
	<DE>		:= "DE   "<description>
	<SQ>		:= "SQ   "<text>

*/

#include <meme.h>

/* size of memory chunks to allocate for sequence data */
#define RCHUNK 100

/* maximum size of sequence description text string */
#define MAXDELEN 10000

/* local types */
typedef enum {FASTA, SWISSPROT} FORMAT_TYPE;

/* local functions */
static SAMPLE *create_sample(
  char *alpha,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence  	/* the sequence */
);

/* local functions */
static long read_sequence_data(
  FILE *data_file,		/* data file of sequences */
  char **sequence, 		/* sequence data */
  char *name			/* name of sequence */
);
static int read_sequence_de(
  FILE *data_file,		/* data file of sequences */
  char **descriptor		/* sequence descriptor */
);

/* local variables */
#define DATA_HASH_SIZE 100000
/* 
   The hash table will be appended to if read_seq_file is called more
   than once.  This is so that positive sequences the negative sequence file
   will be ignored.
*/ 
static HASH_TABLE ht_seq_names = NULL;	/* hash of dataset seq names */

/**********************************************************************/
/*
	read_seq_file

	Open a sequence file and read in the sequences. 
	Returns a dataset->

	setup_hash_alphabet must have been called prior to first call to this.
*/
/**********************************************************************/
extern DATASET *read_seq_file(
  char *file_name,		/* name of file to open */
  char *alpha,			/* alphabet used in sequences */
  BOOLEAN use_comp,		/* use complementary strands, too */
  double seqfrac		/* fraction of input sequences to use */
)
{
  int i, j;
  FILE *data_file;		/* file with samples to read */
  char *sample_name;		/* name of sample read */
  char *sample_de;		/* descriptor text for sample */
  char *sequence;		/* sequence read */
  long length;			/* length of sequence */
  BOOLEAN error=FALSE;		/* none yet */
  SAMPLE *sample;		/* sample created */
  DATASET *dataset;		/* dataset created */
  int n_samples=0;		/* number of samples read */
  double *seq_weights=NULL;	/* sequence weights */
  int n_wgts=0;			/* number of sequence weights given */

  /* create a hash table of sequence names */
  if (!ht_seq_names) ht_seq_names = hash_create(DATA_HASH_SIZE);

  /* create a dataset */
  dataset = (DATASET *) mymalloc(sizeof(DATASET));
  dataset->alength = strlen(alpha);
  dataset->alphabet = alpha;

  /* open data file */ 
  if (file_name == NULL) {
    fprintf(stderr, "You must specify a data file or `stdin'.\n");
    exit(1);
  } else if (strcmp(file_name, "stdin")) {
    data_file = fopen(file_name, "r"); 
    if (data_file == NULL) {
      fprintf(stderr, "Cannot open file `%s'.\n", file_name);
      exit(1);
    }
  } else {
    data_file = stdin;
  }

  /* initialize maximum length of sequences */
  dataset->max_slength = 0;
  dataset->min_slength = 10000000;

  dataset->n_samples = 0;	/* no samples yet */
  dataset->samples = NULL;	/* no samples */

  while (read_sequence(data_file, &sample_name, &sample_de, &sequence, 
    &length)) {

    /* skip sequence if an error occurred */
    if (length < 0) continue;

    /* parse weights if given; make (more than enough) room in array */
    if (strcmp(sample_name, "WEIGHTS")==0) {
      double wgt; 
      char *wgt_str = sample_de;
      Resize(seq_weights, n_wgts+(int)strlen(wgt_str), double);
      while (sscanf(wgt_str, "%lf", &wgt) == 1) {
        if (wgt <= 0 || wgt > 1) {
	  fprintf(stderr, 
            "Weights must be larger than zero and no greater than 1.\n");
	  exit(1);
        }
        seq_weights[n_wgts++] = wgt;			/* save weight */
        wgt_str += strspn(wgt_str, "      ");		/* skip white */
        wgt_str += strcspn(wgt_str, "     ");		/* skip token */
      }
      myfree(sample_name);
      myfree(sample_de);
      myfree(sequence);
      continue;
    }

    /* ignore duplicate (same sample name) sequences */ 
    if (hash_lookup(sample_name, 0, ht_seq_names)) {
      fprintf(stderr, "Skipping sequence '%s'.\n", sample_name);
      myfree(sample_name);
      myfree(sample_de);
      myfree(sequence);
      continue;
    }
    hash_insert(sample_name, 0, ht_seq_names);  /* put name in hash table */

    n_samples++;

    /* see if sequence will be used in random sample; store it if yes */
    if (drand48() >= 1 - seqfrac) {

      /* create the sample */
      sample = create_sample(alpha, length, sample_name, sequence);
      if (sample == NULL) {error = TRUE; continue;}

      /* record maximum length of actual sequences */
      dataset->max_slength = MAX(sample->length, dataset->max_slength);
      dataset->min_slength = MIN(sample->length, dataset->min_slength);

      /* put the sample in the array of samples */
      if ((dataset->n_samples % RCHUNK) == 0) {
        Resize(dataset->samples, dataset->n_samples + RCHUNK, SAMPLE *);
      }
      dataset->samples[dataset->n_samples++] = sample;

    }

  } /* sequences */
  if (length < 0) error = TRUE;			/* read_sequence error */

  /* resize the array of samples */
  if (dataset->n_samples) Resize(dataset->samples, dataset->n_samples, SAMPLE*);

  /* check that datafile contained at least one sample */
  if (!error) {
    if (n_samples == 0) {
      fprintf(stderr, "No sequences found in file `%s'.  Check file format.\n",
	file_name);
      error = TRUE;
    } else if (dataset->n_samples == 0) {
      fprintf(stderr, 
        "No sequences sampled.  Use different seed or higher seqfrac.\n");
      error = TRUE;
    }
  }

  /* exit if there was an error */
  if (error) exit(1);

  /* calculate the prior residue frequencies and entropies 
     and |D|, size of the dataset */
  /* tlb; 5/9/97 wgt_total_res and weighted res_freq */
  dataset->res_freq = NULL;
  Resize(dataset->res_freq, dataset->alength, double);
  for (i=0; i<dataset->alength; i++) { dataset->res_freq[i] = 0; }
  dataset->total_res = 0;
  dataset->wgt_total_res = 0;
  for (i=0; i<dataset->n_samples; i++) {		/* sequence */
    long slen = dataset->samples[i]->length;
    double sw = dataset->samples[i]->sw = (n_wgts > i ? seq_weights[i] : 1);
    dataset->total_res += slen;
    dataset->wgt_total_res += slen*sw;
    for (j=0; j<dataset->alength; j++) {
      if (use_comp) { 	/* using complementary strand as well */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j]/2.0;
        dataset->res_freq[hash(comp_dna(unhash(j)))] += 
          sw * dataset->samples[i]->counts[j]/2.0;
      } else {		/* not using complementary strands */
        dataset->res_freq[j] += sw * dataset->samples[i]->counts[j];
      }
    }
  }

  /* convert counts to frequencies */
  for (i=0; i<dataset->alength; i++)  
    dataset->res_freq[i] /= dataset->wgt_total_res;

  return dataset;
} /* read_seq_file */

/**********************************************************************/
/*
	create_sample

	Create a sample.

	Returns the sample or NULL on error.

*/
/**********************************************************************/
static SAMPLE *create_sample(
  char *alpha,		/* the alphabet */
  long length,		/* length of sequence */
  char *name,		/* name of sample */
  char *sequence  	/* the sequence */
)
{
  long i, j; 
  SAMPLE *new1;
  int alength = strlen(alpha);			/* length of alphabet */

  /* disallow zero length samples */
  if (length == 0) {
    fprintf(stderr, "\nZero length sequences not allowed. (Sequence `%s').\n",
      name);
    return NULL;
  }

  /* create the record to hold the sample and its associated data */
  new1 = (SAMPLE *) mymalloc(sizeof(SAMPLE));

  /* assign the name and sequence data */
  new1->sample_name = name;
  new1->seq = strdup(sequence);

  /* set up encoded version of sequence and weights */
  new1->res = (char *) mymalloc(length * (int) sizeof(char));
  new1->resic = (char *) mymalloc(length * (int) sizeof(char));
  new1->pYic = (char *) mymalloc(length * (int) sizeof(char));
  new1->weights = (double *) mymalloc(length * (int) sizeof(double));
  new1->not_o = (double *) mymalloc(length * (int) sizeof(double));
  new1->log_not_o = (int *) mymalloc(length * (int) sizeof(int));
  new1->logcumback = (double *) mymalloc((length+1) * (int) sizeof(double));
  new1->nsites = 0;
  new1->sites = NULL;
  new1->z = (double *) mymalloc(length * (int) sizeof(double));
  create_2array(new1->sz, double, 2, length);
  for (i=0; i<length; i++) { 
    int c = (int) sequence[i];
    int e = hash(c);
    new1->res[i] = e;
    new1->resic[length-i-1] = hash(comp_dna(c));
    new1->weights[i] = 1.0;
  }

  /* set up arrays to hold posterior probabilities */
  create_2array(new1->pY, int, 3, length);

  /* set up array to hold character counts) */
  new1->counts = (double *) calloc(alength, (int) sizeof(double));

  /* get character counts */
  for (i=0; i<length; i++) {
    int c = new1->res[i];
    if (c<alength) {				/* normal letter */
      new1->counts[c]++;
    } else {					/* 'X' */
      for (j=0; j<alength; j++) new1->counts[j] += 1.0/alength;
    }
  }

  /* record length of sequence */
  new1->length = length;

  return new1;
} /* create_sample */

/**********************************************************************/
/*
	read_sequence

	Read a single sequence from the data file.
	Returns FALSE on EOF or bad sequence data.

	Supports FASTA and SWISSPROT formats.

	Note:
	setup_hash_alphabet must have been called prior to first call to this.
*/
/**********************************************************************/
extern BOOLEAN read_sequence(
  FILE *data_file,		/* file containing sequences */
  char **sample_name,		/* unique identifier of sequence */
  char **sample_de,		/* descriptor of sequence */
  char **sequence,		/* protein or DNA letters */
  long *length			/* length of sequence */
)
{
  int i, c;
  int msn = 2 * MSN;			/* define maximum internal name size */
  FORMAT_TYPE format = FASTA;		/* assume FASTA format */
  char *name = NULL;

  /* skip anything until first sample name */
  c = ' '; 
  while(c != EOF) { 
    if((c=fgetc(data_file)) == '>') {	/* FASTA format */
      break;
    } else if (c == 'I') {  		/* swiss-prot format? */
      if ((c = fgetc(data_file)) == 'D') {
	format = SWISSPROT;
        c = fgetc(data_file);
        Skip_whi(c, data_file);		/* skip whitespace to name */
	break;
      }
    } 
    Skip_eol(c, data_file);		/* go to end of line */
  }
  if (c==EOF) return FALSE;		/* no more sequences */

  /* get the sample name; truncate to msn characters */
  Resize(name, msn+1, char);
  /* read to first blank/tab/ or end of line/file */
  for (i=0; (c=fgetc(data_file))!=EOF; ) {
    if ((c==' ') && i==0) {	/* skip blanks until name starts */
      continue;
    } else if (c==' ' || c=='\t' || c=='\n' || c=='\r') {
      break;				/* blank or nl ends name */
    } else if (i < msn) {
      name[i++] = c;			/* non-blank: add to name */
    }
  }
  name[i] = '\0';

  /* read in description */
  *sample_de = NULL;			/* in case no description found */
  if (format == FASTA) {
    if (c != '\n' && c != '\r' && c != EOF) { 
       Skip_whi(c, data_file);		/* skip whitespace to name */
       (void) read_sequence_de(data_file, sample_de);
    }
  } else if (format == SWISSPROT) {
    Skip_eol(c, data_file);				/* go to end of line */
    while (c != EOF) {					/* read all DE lines */
      if ((c=fgetc(data_file)) == 'D') {	 	/* start of DE line? */
	if ((c=fgetc(data_file)) == 'E') {
	  c=fgetc(data_file);
          Skip_whi(c, data_file);			/* skip white space */
	  (void) read_sequence_de(data_file, sample_de);/* read description */
          c = '\n';					/* at end of line */
        }
      } else if (c == 'S') { 				/* start of SQ line? */
	if ((c=fgetc(data_file)) == 'Q') {
          Skip_eol(c, data_file);			/* go to end of line */
          break;
        }
      } else {
        Skip_eol(c, data_file);				/* go to end of line */
      }
    } /* read all DE lines */
  } /* format */

  /* read in the actual sequence data */
  *length = read_sequence_data(data_file, sequence, name);

  /* sequence had bad data */
  if (*length < 0) {
    myfree(name);
    myfree(*sample_de);
    return FALSE;
  }

  *sample_name = name;
  /* insure that the description string exists */
  if (*sample_de == NULL) {
    Resize(*sample_de, 1, char);
    (*sample_de)[0] = '\0';
  }

  return TRUE;
} /* read_sequence */

/**********************************************************************/
/*
	read_sequence_data

	Read the sequence data into a dynamic array.
	Converts to upper case.
	Checks for illegal sequence characters.

	Returns the length of the sequence or -1 on error.
*/
/**********************************************************************/
static long read_sequence_data(
  FILE *data_file,		/* data file of sequences */
  char **sequence, 		/* sequence data */
  char *name			/* name of sequence */
)
{
  long length; 
  int c;
  char *seq = NULL;

  /* 
    read sample sequence 
  */
  for(length=0; (c=fgetc(data_file))!=EOF; ) {
    if (c == '>') { 			/* end of FASTA sequence */
      ungetc(c,data_file); 
      break; 
    } else if (c == '/') {		/* end of SWISSPROT sequence */
      c = fgetc(data_file);
      if (c != '/') {
        fprintf(stderr, "\nError reading SWISSPROT database.\n");
        exit(1);
      }
      while ((c=fgetc(data_file))!=EOF && c != '\n') ;	/* find end of line */
      break;
    } else {
      if (isspace(c)) continue;			/* remove whitespace from seq */
      if(islower(c)) c = toupper(c);		/* convert to uppercase */
      if (hash(c) == -1) {			/* illegal character? */
        if (c=='*' || c=='.') continue;		/* remove *'s and .'s for MEME*/
        fprintf(stderr, "\nIllegal character `%c' in sequence %s.\n", c, name);
        fprintf(stderr, "Change alphabet or fix data file.\n");
        if (seq != NULL) myfree(seq);
        return -1;
      }
      /* alocate space for ASCII version of sequence */
      if ((length % RCHUNK) == 0) {
        Resize(seq, length+RCHUNK, char);
      }
      seq[length++] = c;
    } /* read sample sequence */
  }
  /* put NULL at the end of the sequence */
  Resize(seq, length+1, char);
  seq[length] = '\0';

  *sequence = seq;

  return length;
} /* read_sequence_data */

/**********************************************************************/
/*
	read_sequence_de

	Read the sequence descriptor into a dynamic array.
*/
/**********************************************************************/
static int read_sequence_de(
  FILE *data_file,		/* data file of sequences */
  char **descriptor		/* sequence descriptor; if not NULL,
				   string will be appended */
)
{
  int length, c;
  char *de = *descriptor;

  /* find current length of descriptor */
  length = 0;
  if (de != NULL) {			/* continuation of prev. descriptor */
    length = strlen(de);		/* old length */
    /* alocate space for descriptor string */
    Resize(de, (1+length/RCHUNK)*RCHUNK, char);
    de[length++] = ' ';			/* put a blank at the end */
  }

  /* read in sample descriptor */
  while ( (c=fgetc(data_file)) != EOF) {
    if (c == '\n') { 
      break; 
    } else if (length <= MAXDELEN) {
      /* alocate space for descriptor string */
      if ((length % RCHUNK) == 0) Resize(de, length+RCHUNK, char);
      de[length++] = c;			/* append character to de */
    }
  }
  /* put a NULL at the end of the de */
  Resize(de, length+1, char);
  de[length] = '\0';

  *descriptor = de;

  return length;
}
