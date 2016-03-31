/*
 * $Id: meme.c,v 1.2.4.1 2006/01/25 08:06:03 tbailey Exp $
 * 
 * $Log: meme.c,v $
 * Revision 1.2.4.1  2006/01/25 08:06:03  tbailey
 * Print "CPU: " even if UNIX not defined so output file will pass automatic
 * tests.
 *
 * Revision 1.2  2005/10/25 19:06:39  nadya
 * rm old macro for Header, all info is taken care of by Id and Log.
 *
 * Revision 1.1.1.1  2005/07/29 17:19:58  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/************************************************************************
*								       	*
*	MEME							       	*
*	Author: Timothy L. Bailey				       	*
*									*
*	Copyright							*
*	(1994 - 2000) The Regents of the University of California.	*
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
*								       	*
************************************************************************/

/* 5-26-00 tlb; initialize model->cons0 in init_model */
/* 10-3-99 tlb; replace nrdb with back */
/* 7-29-99 tlb; set model->pal correctly in init_model */
/* 7-02-99 tlb; in erase normalize z_ij to 1.0 so erasing will be complete */
/* 7-02-99 tlb; in erase run just e_step of em */
/* 6-28-99 tlb; remove sigma from model */
/* 8-7-97 tlb; removed Mcm stuff from erase */

/*	
	meme <datafile> [options]

	Reads in a sequence file (Pearson/FASTA format).
	Runs EM algorithm with selected values of W and lambda.

	<datafile>	set of samples: [>id sequence]+
*/

	
#define DEFINE_GLOBALS 
#include "meme.h"

/* FIXME ??? */
#ifndef PARALLEL
#define mpMyID() 0
#endif

/* external subroutines */
extern double sqrt(double x);

/* local subroutines */
static void erase(
  DATASET *dataset,			/* the dataset */
  MODEL *model	 			/* the model */
);
/* added by M.H. */
static void readSecondaryStructurePrior(
  DATASET *dataset,			/* the dataset */
  MODEL *model,	 			/* the model */
  char *secondaryStructureFilename,	/* input filename */
  double secondaryStructurePseudocount  /* pseudocount for prior prob */
);

static BOOLEAN save_candidate(
  MODEL *model,				/* final model */
  DATASET *dataset,			/* the dataset */
  S_POINT *s_point,			/* starting point */
  CANDIDATE *candidates,		/* candidate for best model of width */
  double best_sig			/* best significance */
);
static BOOLEAN init_model(
  S_POINT *s_point,			/* the starting point */
  MODEL *model,				/* the model to intialize */
  DATASET *dataset,			/* the dataset */
  int imotif				/* motif number */
);

BOOLEAN no_print = FALSE;	/* turn off printing if parallel and not main */

/**********************************************************************/
/*
	main 
*/
/**********************************************************************/

extern int main(
  int argc,
  char *argv[]
)
{
  int i, imotif;
  DATASET *dataset;		/* the dataset */
  DATASET *neg_dataset;		/* the dataset of negative examples */
  MODEL *model;			/* the model */
  MODEL *neg_model;		/* the model of negatives */
  MODEL *best_model;		/* the best model found (so far) */
  MODEL *scratch_model;		/* a scratch model */
  CANDIDATE *candidates = NULL;	/* list of candidate models */
  int n_starts = 0;		/* number of starting points */
  int nmotifs;			/* number of motifs to find */
  double stop_time=0, m1time=0;	/* time when stopped, for motif 1 (secs) */

#ifdef PARALLEL
  int start_start, incr;
  /* Initialize MPI. */
  mpInit(&argc, &argv);
  /* turn off printing if parallel and not the main processor */
  no_print = (mpMyID() != 0);
  incr = mpNodes();
  /*fprintf(stderr, "Running on %d nodes...\n", incr);*/
#endif /* PARALLEL */

#ifdef debug_ieee
  ieee_handler_setup("common");
#endif 

  (void) myclock();		/* record CPU time */

  /* print the command line */
  if (VERBOSE) {
    argv[0] = "meme";
    for (i=0; i<argc; i++) printf("%s ", argv[i]);
    printf("\n\n");
    fflush(stdout);
  }
#ifdef UNIX
  if (!no_print && VERBOSE) system("echo ''; echo CPU: `hostname`; echo '' ");
#else
  if (!no_print && VERBOSE) printf("\nCPU: unknown\n\n");
#endif /* UNIX */

  /* initialize everything from the command line */
  init_meme(argc, argv, &model, &best_model, &scratch_model, &neg_model, 
    &dataset, &neg_dataset);
  nmotifs = dataset->nmotifs;		/* number of motifs to find */
  Resize(candidates, model->max_w+2, CANDIDATE);        /* motif candidates; max_w+2 in case -pal given */

  /* describe the dataset */
  print_dataset_summary (dataset);

  /* describe the command */
  print_command_summary(model, dataset);
 
  /* print explanation of how to interpret the results */
  if (DOC) print_meme_doc(dataset->meme_directory);

  if (!NO_STATUS) {
    fprintf(stderr, "\nseqs=%6d, min=%4ld, max=%5ld, total=%9d\n",
      dataset->n_samples, dataset->min_slength, dataset->max_slength,
      dataset->total_res);
  }

	/* added by M.H. */
   /* if a filename with sec struct information is given this function computes a new prior for P(Z_ij = 1)
		that is not uniformly distributed; 
		if no file is given, standard MEME is used */
   if (dataset->secondaryStructureFilename != NULL) {
		printf("use secondary structure information from: %s\n", dataset->secondaryStructureFilename); fflush(stdout);
      readSecondaryStructurePrior(dataset, model, dataset->secondaryStructureFilename, dataset->secondaryStructurePseudocount);
	 	printf("done\n"); fflush(stdout);
	}	


  /*  Find a concept and erase it loop */
  for (imotif=1; imotif<=nmotifs; imotif++) {
    S_POINT *s_points = NULL;		/* array of starting points */
    int i_start;			/* index in s_points */
    char *e_cons = dataset->p_point->e_cons0[imotif-1];	/* consensus sequence */
    int best_w = 0;			/* best width */
    double best_sig;			/* motif significance */
    double w_factor = sqrt(2.0);	/* factor separating sampled widths */
    int iter = 0;                       /* total number of EM iterations */

    if (!NO_STATUS) fprintf(stderr, "\nmotif=%d\n", imotif);

    /* known motif has been given */
    if (dataset->nkmotifs > 0) {
      model->min_w = model->max_w = dataset->motifs[imotif-1].width;
      model->min_nsites = model->max_nsites = dataset->motifs[imotif-1].pos;
    }

    /* set up the array of starting points for EM */
    s_points = get_starts(dataset, model, e_cons, w_factor, &n_starts);

    /* leave loop if no starts found */
    if (n_starts == 0) break;

    /* tag each candidate width as unused; max_w+1 in case -pal given */
    for (i=0; i<=model->max_w+1; i++) {
      candidates[i].sig = BIG;
      candidates[i].s_point = NULL;
    }

    /* run EM on each start and save best final model for each final width */
    best_sig = BIG;				/* best motif significance */

#ifdef PARALLEL
    /* Check whether to parallelize this loop. */
    start_start = mpMyID();
    /* Make sure everybody has something to do. */
    if (start_start >= n_starts) start_start = n_starts-1;
    /* Divide the various starting points among processors. */
    for (i_start=start_start; i_start < n_starts; i_start += incr) {
#else
    for (i_start=0; i_start<n_starts; i_start++) {
#endif /* PARALLEL */
      S_POINT *s_point = s_points+i_start;	/* current starting point */
#ifdef DEBUG
      double s_point_time = myclock()/1E6;
#endif /* DEBUG */

      /* initialize the model from the starting point */
      if (! init_model(s_point, model, dataset, imotif)) continue;

      /* Count iters per loop. */
      model->iter = 0;

      /* Run EM from the starting model */
      em(model, dataset);

      /* Keep track of the total number of EM iterations. */
      iter += model->iter;

      /* store model as a candidate for best model;
	 save model if best so far */
      if (save_candidate(model, dataset, s_point, candidates, best_sig)) {
	best_w = model->w;			/* new best width */
	swap(model, best_model, MODEL *);	/* save model as best_model */
	best_sig = candidates[best_w].sig;	/* new best significance */
      }

    } /* starting point loop */ 

    /* found the best model */
    swap(model, best_model, MODEL *);

#ifdef PARALLEL
    /* Copy the consensus sequence into the model. */
    store_consensus(model, candidates);

    /* Do the reduction. */
    reduce_across_models(model, dataset->alength);
    /*fprintf(stderr, "%d: past reduction\n", mpMyID()); fflush(stderr);*/
#endif /* PARALLEL */

    /* quit if model has too few sites */
    if (model->nsites_dis < MINSITES) break;
    /*fprintf(stderr, "%d: past few sites\n", mpMyID()); fflush(stderr);*/

    /* quit if model fails E-value test */
    if (model->logev > log(dataset->evt)) break;
    /*fprintf(stderr, "%d: past E-value\n", mpMyID()); fflush(stderr);*/

    /* Store the total number of EM's in the model. */
    model->iter = iter;

    /* ERASE the site and starts */
    /*fprintf(stderr, "%d: at erase\n", mpMyID()); fflush(stderr);*/
    erase(dataset, model);
    /*fprintf(stderr, "%d: past erase\n", mpMyID()); fflush(stderr);*/

    /* calculate negative model by doing one EM iteration on negative ex's */
    if (neg_model) {
      /* copy motif model to negative model */
      copy_model(model, neg_model, dataset->alength);
      /* get negative model */
      neg_dataset->maxiter = 1;			/* 1 iteration of em */
      em(neg_model, neg_dataset);
      /* ERASE the site and starts */
      erase(neg_dataset, neg_model);
    }

    /* print results */
    /*fprintf(stderr, "%d: at print results\n", mpMyID()); fflush(stderr);*/
    if (!no_print) print_results(dataset, neg_dataset, model, neg_model,
      candidates);
    /*fprintf(stderr, "%d: past print results\n", mpMyID()); fflush(stderr);*/

    /* stop if out of time */
    if (dataset->max_time && imotif<nmotifs) { 	/* considering another motif */
      stop_time = myclock()/1E6;		/* current time */
#ifdef PARALLEL
      /* use stop_time from process 0 to avoid it stopping while others continue*/
      /*fprintf(stderr, "%d: time broadcast size %d\n", mpMyID(),sizeof(stop_time)); fflush(stderr);
	  */
      mpBroadcast((void *)&stop_time, sizeof(stop_time), 0);
      /*fprintf(stderr, "%d: past time broadcast\n", mpMyID()); fflush(stderr); */
#endif
      /* record time if this is motif 1 */
      if (imotif == 1) m1time = stop_time;
      if ((dataset->max_time - stop_time) < m1time) {
        ++imotif;				/* this motif OK */
        break;
      }
    } /* check time */
    
    /*fprintf(stderr, "%d: past check time\n", mpMyID()); fflush(stderr); */
  } /* nmotifs loop */ 
  --imotif;				/* number of motifs found */

  /* print the motif block diagrams using all the motifs */
  if (!no_print && imotif) 
     print_summary(model, dataset, los, imotif, pv, stdout);

  /* print reason for stopping */
  printf("\n"); PSTARS;
  if (n_starts == 0) {
    printf("Stopped because couldn't find any more starting points for EM.\n");
  } else if (imotif == nmotifs) {
    printf("Stopped because nmotifs = %d reached.\n", nmotifs);
  } else if (dataset->max_time && (dataset->max_time - stop_time) < m1time) {
    printf("Stopped because would probably run out of time (%.2f secs).\n", 
      dataset->max_time);
  } else if (model->nsites_dis < MINSITES) {
    printf("Stopped because next motif has fewer than %d sites.\n",
      MINSITES);
  } else {
    printf("Stopped because motif E-value > %8.2e.\n", dataset->evt);
  }
  PSTARS;
  fflush(stdout);

#ifdef UNIX
  if (!no_print) system("echo ''; echo CPU: `hostname`; echo '' ");
#else
  if (!no_print) printf("\nCPU: unknown\n\n");
#endif /* UNIX */
  PSTARS;

  if (!NO_STATUS) fprintf(stderr, "\n");

#ifdef PARALLEL
   mpFinalize();
#endif

	
	/* added by M.H. */
   if ((dataset->secondaryStructureFilename != NULL) && (model->mtype == Tcm)) {
	  printf("max adjustment of pseudocount: %f\n",MAXADJUST);
	}

  return(0);
} /* main */



/* added by M.H. */
/**********************************************************************/
/*
	readSecondaryStructurePrior

     read the secondary structure information (EF or PU) from the filename 
	  to use these as prior probabilities sigma_ij = P(Z_ij = 1)
	  this function stores the EF/PU values and computes the the sigma array for each sample

*/
/**********************************************************************/
static void readSecondaryStructurePrior(
  DATASET *dataset,			/* the dataset */
  MODEL *model,	 			/* the model */
  char *secondaryStructureFilename,	/* input filename */
  double secondaryStructurePseudocount  /* pseudocount for prior prob */
)
{
  int i, j;
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* the sequences */
  FILE *data_file;				
  int line_count = 0;         			/* number of read lines */		
  int maxLen = dataset->max_slength;	/* for determining buffer length */
  int bufLen = 0;								/* buffer length */

  /* open file */
  data_file = fopen(secondaryStructureFilename, "r"); 
  if (data_file == NULL) {
    fprintf(stderr, "Cannot open file %s\n", secondaryStructureFilename);
    exit(1);
  }

  /* for reading the file - we expect 6 chars for each value + ';' + some space for the sequence name (1000 chars) */	
  bufLen = maxLen*7 + 1000;
  char *line = (char*) malloc( bufLen * sizeof(char));        

  /* read the file and store the values in ss_value for each sample */	
  for (;;) {
    char *word;                        /* for tokenize */
    char *seq_name;                    /* name of the sequence */
    double prob;                       /* for conversion string to double */
    int count = 0;                     /* count number of read values */
 
    fgets(line, bufLen, data_file);
    if (feof(data_file)) 
      break;
    /* remove \n at line end */
	 if (line[strlen(line)-1] == '\n') {	
	 	line[strlen(line)-1] = '\0'; 
	 }else{
	 	fprintf(stderr, "ERROR in readSecondaryStructurePrior: line %d contains more than %d characters\n", line_count+1, bufLen);				
		exit(1);
	 }
	 
	 /* skip empty lines */
	 if (strlen(line) == 0) 
	    continue;
	 
	 /* avoid reading more lines than there are sequences */
	 if (line_count >= n_samples) {
	 	fprintf(stderr, "ERROR in readSecondaryStructurePrior: read line %d but there are only %d sequences\n", line_count+1, n_samples);				
		exit(1);
	 }	 
	 
	 /* tokenize this line */
	 word = strtok(line,";");		/* first token is sequence name */
	 seq_name = word;
	 /* make sure that we read the values for the correct sequence */
	 if (strcmp(samples[line_count]->sample_name, seq_name) != 0) {
      fprintf(stderr, "ERROR in readSecondaryStructurePrior: expect secondary structure information for sequence %s but read for sequence %s (maybe the file is shuffled)\n", samples[line_count]->sample_name, seq_name);
      exit(1);	 
	 }
	 
    int lseq = samples[line_count]->length;     /* seq length */
	 samples[line_count]->ss_value = (double*) malloc (lseq * sizeof(double));  		/* alloc space */
    double *ss_value = samples[line_count]->ss_value;        /* ss_value_ij */

	 /* read remaining tokens (lines that are empty are ignored since they have only one token (NULL)) */
    for (;;) {
      word = strtok(NULL, ";");
	 	if (word == NULL)  
			break;
		prob = atof(word);
	   if (count < lseq) 		/* store current value but avoid writing out of range */
			ss_value[count] = prob;
		count++;						/* count number of read values */

    }
	 if (count != lseq - model->w + 1) {
 		fprintf(stderr, "ERROR in readSecondaryStructurePrior: expect %d secondary structure values for sequence %s but %d values are read (maybe you set the wrong motif length -w)\n", 
			lseq-model->w+1, seq_name, count);
    	exit(1);	 
    }
	  
    line_count	++;
	 
  }
  fclose(data_file);

  /* check if secondary structure values for each sequences are read */
  if (line_count != n_samples) {
    fprintf(stderr, "ERROR in readSecondaryStructurePrior: %d lines read but there are %d sequences\n", line_count, n_samples);				
    exit(1);
  }	 


  /* allocate space for sigma and compute sum and max of ss_value[i] */
  for (i=0; i<n_samples; i++) {
  	 SAMPLE *s = samples[i];
	 int lseq = s->length;

	 s->sigma = (double*) malloc (lseq * sizeof(double));  		/* alloc space */
	 s->intlog_sigma = (int*) malloc (lseq * sizeof(int));  		/* alloc space */

	 /* compute \sum ss_value[i] and maximum ss_value */	
	 s->sum_ss_value = 0;
	 s->max_ss_value = -1;
    for (j=0; j <= lseq - model->w; j++) {
	 	s->sum_ss_value += s->ss_value[j];
		if (s->ss_value[j] > s->max_ss_value) {
			s->max_ss_value = s->ss_value[j];
		}
	 }

	 /* compute sigma, LOG(sigma) and INT_LOG(sigma) 
	 	 INT_LOG(sigma) is used by all three models during the start point search 
		 LOG(sigma) is used by OOPS since the prior remains constant during EM 
		 for ZOOPS sigma remains constant but the LOG have to be taken in each EM iteration since sigma is multiplied by gamma^(t) 
		 for TCM sigma is computed in each iteration because lambda changes
		 
	 */
    double *sigma = s->sigma;					 		 /* sigma_ij */
    int *intlog_sigma = s->intlog_sigma; 		 /* INT_LOG(sigma_ij) */
	 double sum = s->sum_ss_value + ((lseq - model->w + 1) * secondaryStructurePseudocount);  		 /* \sum ss_value[i] + m*pseudocount */

	 /* in case of OOPS precompute the LOG(sigma[i]) since sigma remains constant during EM */
	 if (model->mtype == Oops) {
	    s->log_sigma = (double*) malloc (lseq * sizeof(double));		 /* alloc space */
	 }
    double *log_sigma = s->log_sigma;  		/* sigma_ij */
	 if (VERBOSE) {
		 printf("prior prob for seq %s: ",samples[i]->sample_name);
		 fflush(stdout);
	 }
	 
    for (j=0; j <= lseq - model->w; j++) {
	    sigma[j] = (s->ss_value[j] + secondaryStructurePseudocount) / sum;
	    intlog_sigma[j] = INT_LOG(sigma[j]);
	    if (model->mtype == Oops) {
	   	 log_sigma[j] = LOG(sigma[j]);
	    }

		 if (VERBOSE) 
		    printf("%1.4f;",sigma[j]); 
	 }
	 if (VERBOSE) 
		 printf("\n");
  }	


} /* readSecondaryStructurePrior */






/**********************************************************************/
/*
	erase

        For all models:
	  Reset the weights of the letters to probabilisticaly "erase"
	  the ones which occur in sites already found.

*/
/**********************************************************************/
static void erase(
  DATASET *dataset,			/* the dataset */
  MODEL *model	 			/* the model */
)
{
  int i, j, k;
  int n_samples = dataset->n_samples;		/* number of sequences */
  SAMPLE **samples = dataset->samples;		/* the sequences */
  int w = model->w;				/* width of motif */

  /* 
    Set z from the maxima stored in the learned model.
  */
    /*fprintf(stderr, "%d: at set_z\n", mpMyID()); fflush(stderr); */
  set_z(model, dataset);
    /*fprintf(stderr, "%d: past set_z\n", mpMyID()); fflush(stderr); */

  /* 
     z_ij is taken as the probability of a site occurring at i,j.
     The probability of a position being in a site is taken
     as the maximum of the z_ij for sites containing (overlapping) it.
     w_ij is set to 1-max(z_ij) times its previous value which
     reflects the independence assumption among motifs.
  */
  for (i=0; i<n_samples; i++) 		{	/* sequence */
    double *weights = samples[i]->weights;	/* w_ij */
    int lseq = samples[i]->length;		/* seq length */
    double *zi = samples[i]->z;			/* z_ij */

    /*if (mpMyID()==0)fprintf(stderr, "0: i loop %d lseq %d\n", i, lseq); fflush(stderr); */

    if (lseq < w) continue;			/* sample too short for motif */

    for (j=0; j<lseq; j++) {			/* position */
      double max_z = 0.0;
      /*if (mpMyID()==0)fprintf(stderr, "0: j loop %d\n", j); fflush(stderr); */
      /* find largest probability that site overlaps this position */
      for (k=MAX(0,j-w+1); k<=j && k<lseq-w+1; k++) {
        /*if (mpMyID()==0)fprintf(stderr, "0: k loop %d\n", k); fflush(stderr);
		*/
	max_z = MAX(max_z, zi[k]);
      }
      max_z = MIN(1.0, max_z);			/* fix roundoff errors */
      /* update the probability that position not in a site */
      weights[j] *= 1.0 - max_z;
    }
  }

  if (PRINT_W) print_wij(dataset);
} /* erase */

/**********************************************************************/
/*
	init_model

	Initialize a model from a starting point.

	Returns false if starting point was not valid.
*/
/**********************************************************************/
static BOOLEAN init_model(
  S_POINT *s_point,			/* the starting point */
  MODEL *model,				/* the model to intialize */
  DATASET *dataset,			/* the dataset */
  int imotif				/* motif number */
)
{
  int w0; 
  
  /* skip if no good starting points found for w0, nsites0 */
  if (s_point->score == LITTLE) { return FALSE; }

  /* initialize the new motif */
  strcpy(model->cons0, s_point->cons0);
  w0 = model->w = model->pw = s_point->w0;
  init_theta(model->theta, s_point->e_cons0, w0, dataset->map,dataset->alength);

  /* initialize lambda */
  model->lambda = MIN(s_point->nsites0/wps(dataset, w0), 1);
  model->pal = dataset->pal;

  /* initialize prior estimate of number of sites */
  model->psites = s_point->nsites0;

  if (PRINTALL) {
    printf("component %2d: lambda= %8.6f ps= %8.0f\n", 
      1, model->lambda, wps(dataset, w0));
    print_theta(0, 2, model->nsites_dis, model->theta, model->w, 0, "", 
      dataset, stdout);
  }

  /* init motif number */
  model->imotif = imotif;
  model->iseq = s_point->iseq;
  model->ioff = s_point->ioff;

  return TRUE;
} /* init_model */

/**********************************************************************/
/*
	save_candidate

	Save the starting point and part of model if it is
	most significant model of its width so far:
		model->sig is smallest

	Returns true if the model is the best so far among all widths.
*/ 
/**********************************************************************/
static BOOLEAN save_candidate(
  MODEL *model,				/* final model */
  DATASET *dataset,			/* the dataset */
  S_POINT *s_point,			/* starting point */
  CANDIDATE *candidates,		/* candidate for best model of width */
  double best_sig			/* best motif significance */
)
{
  int w = model->w;			/* final motif w */
  /* objective function value */
  double sig = dataset->objfun==Pv ? model->logpv : model->logev;

  /* print the results for this w0, nsites0 and THETA */ 
  if (PRINT_STARTS) { 
    printf("\n(start) %3d %6.1f %.*s --> %s ", 
      s_point->w0, s_point->nsites0, s_point->w0, s_point->cons0, model->cons);
    printf("w %3d nsites %4d sig %20.10g\n\n", w, model->nsites_dis, exp(sig));
    fflush(stdout);
  } 

  /* save the results if best so far for this width */ 
  if (sig < candidates[w].sig) {
    candidates[w].s_point = s_point;
    candidates[w].w = w;
    candidates[w].pal = model->pal;
    candidates[w].invcomp = model->invcomp;
    candidates[w].lambda = model->lambda;
    strcpy(candidates[w].cons, model->cons);
    candidates[w].rel = model->rel;
    candidates[w].ll = model->ll;
    candidates[w].sig = sig;
  }

  return (sig < best_sig);
} /* save candidates */

