/*
 * $Id: meme-server.c,v 1.6 2006/01/03 06:37:00 tbailey Exp $
 * 
 * $Log: meme-server.c,v $
 * Revision 1.6  2006/01/03 06:37:00  tbailey
 * Fix "MEME_BINbin" bug in mast-client.txt.
 * Fix indentation in meme-server.c and meme-client.c.
 *
 * Revision 1.5  2005/10/02 00:18:25  nadya
 * revert to mast from mast.csh
 *
 * Revision 1.4  2005/09/13 21:11:36  nadya
 * rm exit line from mast/meme submit scripts
 *
 * Revision 1.3  2005/08/30 22:46:21  nadya
 * update calls to mailer
 *
 * Revision 1.2  2005/08/24 04:37:07  nadya
 * change script layout. Don't assume a given directory,
 * just use full path where needed for files and binaries.
 *
 * Revision 1.1.1.1  2005/07/29 17:17:54  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*
	meme-server

	Usage: meme-server <socket> <command> <qcmd>
		<socket>  socket number to listen on
		<command> command to execute; quote if spaces
		<qcmd>   queueing command including switches and %s 
			  which will be replaced by the qfile name

	If the file sent to the socket contains just the word "ping",
	the program exits with success.

	2/10/96; created; bgrundy and tbailey
*/

/*
Date: Thu, 6 Jul 1995 09:42:01 +0059 (PDT)
From: Allan Snavely <allans@SDSC.EDU>
Subject: Re: sockets at SDSC <Pine.3.89.9506281510.b717-0100000> 
To: Bill Grundy <bgrundy@cs.ucsd.edu>


	Howdy.  I am enclosing two programs supplied by Siamak.
 They realize a simple client server between a Paragon and a 
 workstation.  They don't actually do much; the server (running on
 the // machine) polls a port.  When it gets a message on that port
 from the client it forks a process.  I just tested this between
 ernie (my machine) and xray (both outside the firewall) and it
 worked.

 Here's the server.  Compile icc <file.c> -lm
                     Run a.out -pn open -sz 1  (on the xray)

-------------------------  server follows -----------------------------

Format for input file read from the socket:
----------------------
ADDRESS <return address>
PROGRAM meme
DESCRIPTION <1 line of text>
SWITCHES <switches>
LOGINFO <information to log>
BEGIN
<sequences>


Error messages:
---------------
  10 : Error opening file
  11 : Error reading from socket
  12 : Error writing to file
  13 : Error sending ack to client
  14 : Error binding socket
  15 : Error listening to socket
  16 : Error creating socket
  20 : Header field not found
  21 : Invalid field name
  30 : system error
*/

/***********************************************************************
 * Include files
 ***********************************************************************/
#include "meme-cs.h"

/***********************************************************************
 * Constants
 ***********************************************************************/
#define MAXHEADERLINE 256

/* ParaMEME parameters */
#define MINW 12
#define MAXW 55
#define MAXITER 20
#define TIME "2:00:00"

/*
 system dependent things 
*/

/* returned by system() */
#ifdef sunsparc
  int SYSTEM_OK = -1;
#else
  int SYSTEM_OK = 0;
#endif

/***********************************************************************
 * Type definitions
 ***********************************************************************/
typedef struct HEADER {
  char address[MAXHEADERLINE];
  char program[MAXHEADERLINE];
  char remotename[MAXHEADERLINE];
  char description[MAXHEADERLINE];
  char switches[MAXHEADERLINE];
  char loginfo[MAXHEADERLINE];
} HEADER;


/***********************************************************************
 * Globals 
 ***********************************************************************/

int SOCKET_NUMBER;		/* socket to listen on */
char *MEME;			/* name of meme executable */
char *QCMD;			/* the queueing command */


/***********************************************************************
 *
 * void strlower
 *
 * Downcase a given string.
 *
 ***********************************************************************/
void strlower(char *a_string)
{
  int index;

  for (index = 0; a_string[index] != '\0'; index++)
    a_string[index] = tolower((int)a_string[index]);
}

/***********************************************************************
 *
 * void read_header_item
 *
 * Read one header line from a socket into a HEADER data structure.
 *
 ***********************************************************************/
void read_header_item(HEADER *a_header, int sock, char *fieldname, int field)
{
  int i;
  char buffer[MAXHEADERLINE];
  char *info;
  
  /* Read the line into a local buffer. */
  if (readline(sock, buffer, MAXHEADERLINE) == 0) {
    fprintf(stderr, "job: %d Error reading from socket.\n", (int) getpid());
    fflush(stderr);
    exit(11);
  }
  /*printf("%s (%d)\n", buffer, strlen(buffer));*/

  /* Check to see if this is a ping */
  if (strcmp(buffer, "ping") == 0) {
    exit(0);				/* ping; exit with OK status */
  }

  /* print job number if first field */
  if (field == 1) {
    printf("\njob: %d\n", (int) getpid());
    fflush(stdout);
  }

  /* Search for the field name in the line. */
  if ((info = (char *)strstr(buffer, fieldname)) == NULL) {
    fprintf(stderr, "Field %s not found in %s.\n", fieldname, buffer);
    fflush(stderr);
    exit(20);
  }

  /* Move past the end of the field name. */
  info += strlen(fieldname) + 1;

  /* Get the first character of the field name. */
  /* (Happily, they are mutually exclusive.) */
  switch (fieldname[0]) {
    case 'A' :  /* ADDRESS */
      strcpy(a_header->address, info);
      break;
    case 'P' :  /* PROGRAM */
      strcpy(a_header->program, info);
      break;
    case 'R' :  /* REMOTENAME */
      strcpy(a_header->remotename, info);
      break;
    case 'D' :  /* DESCRIPTION */
      strcpy(a_header->description, info);
      break;
    case 'S' :  /* SWITCHES */
      /* copy removing any nasty semi-colons to prevent security breach */
      for (i=0; info[i]; i++) {
        if (info[i] != ';') {
          a_header->switches[i] = info[i];
        } else {
          a_header->switches[i] = ' ';
        }
      }
      a_header->switches[i] = '\0';
      break;
    case 'L' :  /* LOGINFO */
      strcpy(a_header->loginfo, info);
      break;
    case 'B' :  /* BEGIN */
      /* Don't do anything. */
      break;
    default:
      fprintf(stderr, "Invalid fieldname: %s\n", fieldname);
      fflush(stderr);
      exit(21);
  }
}


/***********************************************************************
 *
 * write_data_file
 *
 * Read in some sequences from a socket and write them to a file.
 *
 ***********************************************************************/
void write_data_file(int sock, char *filename)
{
  int filedescriptor, nread;
  char achar;

  /* Create the data file and open it in write-only mode. */
  if ((filedescriptor = open(filename, 
			     O_WRONLY | O_CREAT | O_TRUNC, 400)) == -1) {
    fprintf(stderr, "Error creating file %s\n", filename);
    fflush(stderr);
    exit(10);
  }
      
  /* Read from the socket until no more bytes come through. */
  for (nread = readn(sock, &achar, 1); nread == 1;
       nread = readn(sock, &achar, 1)) {
    
    /*printf("%c", achar);*/

    /* Echo incoming data to the file. */
    if (write(filedescriptor, &achar, 1) != 1) {
      fprintf(stderr, "Error writing to file %s.\n", filename);
      fflush(stderr);
      exit(12);
    }
  }

  if (nread != 0) {
    fprintf(stderr, "Error reading from socket.\n");
    fflush(stderr);
    exit(11);
  }
     
  /* Close the file. */
  close(filedescriptor);
}

/***********************************************************************
 *
 * receive_file
 *
 * Read in a file from the socket.
 *
 ***********************************************************************/
void receive_file(int sock, char *filename, HEADER *a_header)
{
  /* Read in the header. */
  read_header_item(a_header, sock, "ADDRESS", 1);
  read_header_item(a_header, sock, "PROGRAM", 2);
  read_header_item(a_header, sock, "REMOTENAME", 3);
  read_header_item(a_header, sock, "DESCRIPTION", 4);
  read_header_item(a_header, sock, "SWITCHES", 5);
  read_header_item(a_header, sock, "LOGINFO", 6);
  read_header_item(a_header, sock, "BEGIN", 7);

  /* Echo the sequence data to a file. */
  write_data_file(sock, filename);
}

/***********************************************************************
 *
 * send_ack
 *
 * Send a single-byte ack to the client.
 *
 ***********************************************************************/
void send_ack(int sock)
{
  char achar = '1';

  if (writen(sock, &achar, 1) != 1) {
    fprintf(stderr, "Error sending ack to client.\n");
    fflush(stderr);
    exit(13);
  }
}

/***********************************************************************
 *
 * char *newfilename
 *
 * Creates a new filename, using the process ID for uniqueness.
 *      meme.<pid>.<filename>
 *
 ***********************************************************************/
char *newfilename(char *filename, char * dir)
{
  int length = 6 + 20 + strlen(filename) + 1 + strlen(dir);
  char *name = (char *) mymalloc(length);
 
  /* $dir is a palce holder, it will be filled in a shells cript */
  sprintf(name, "%s/meme.%ld.%s", dir, (long) getpid(), filename);
  return name;
}

/***********************************************************************
 *
 * void make_q_script
 *
 * Write out an NQS submission file for the current job.
 *
 ***********************************************************************/
void make_q_script(HEADER *a_header, char *meme_input, char *qfilename, 
                   char *meme_bin, char *meme_logs)
{
  FILE *qfile;
  char *resfilename, *datefilename;
  long job = (long) getpid(); 		/* get the job number */
  char dcmd[256];			/* date command */
  /*char *dfmt = "-u '+\%d/\%m/\%y \%H:\%M:\%S'";*/	/* date format */
  char *dfmt = "-u '+%d/%m/%y %H:%M:%S'";	/* date format */
  char *type, *text, *npstr;
  int nprocs = 0;			/* assume non-parallel */
  char *p0;

  
  #include <errno.h>

  /* create temporary file names */
  resfilename = newfilename("results", meme_logs);
  datefilename = newfilename("date", meme_logs);

  /* get the type of output: html or ascii */
  type = strstr(a_header->description, "(Use web browser to view results)") ?
    "-html" : "";

  /* send MAST as text file if MEME is */
  text = strstr(a_header->description, "(Use web browser to view results)") ? 
    "" : "-text";

  /* get number of processes to use if parallel MEME */
  if ((npstr = strstr(MEME, "-p ")) != NULL) { 		/* nprocs given */
    sscanf(npstr, "-p %d", &nprocs);
  }

  /* create a script to:
	1) execute MEME
	2) mail the result of MEME
	3) execute and mail results of MAST
  */

  /* Open the script file. */
  if ((qfile = fopen(qfilename, "w")) == NULL) {
    fprintf(stderr, "Error opening file %s.\n", qfilename);
    fflush(stderr);
    exit(10);
  }

  /* script header */
  fprintf(qfile, "#!/bin/csh\n");
  fprintf(qfile, "# This script was automatically generated ");
  fprintf(qfile, "by the MEME server (meme-server.c).\n\n");


  /* change to LOGS directory (needed by qsub) */
  fprintf(qfile, "set dir = %s\n", meme_logs);
  fprintf(qfile, "set bin = %s\n", meme_bin);

  /* save the submission time */
  sprintf(dcmd, "date %s > %s\n", dfmt, datefilename);
  system(dcmd);					/* get submission time */

  /* save the start time */
  fprintf(qfile, "set t1 = `date %s`\n", dfmt);

  /* run MEME and save the results in a file */
  if ((p0 = strstr(MEME, "-p ")) != NULL) {	/* nprocs given */
    int flen;
    p0 += 3;
    flen = (int)(p0 - MEME);
    fprintf(qfile, "%*.*s\"%s\"\\\n", flen, flen, MEME, p0);
  } else if (strstr(MEME, "-p") != NULL) { 
    fprintf(qfile, "%s %d\\\n", MEME, nprocs);
  } else {
    fprintf(qfile, "%s\\\n", MEME);
  }
  fprintf(qfile, "  %s \\\n", meme_input);
  fprintf(qfile, "  %s \\\n", a_header->switches);      /* overridable */
  fprintf(qfile, "  -nostatus -maxiter %d \\\n", MAXITER);
  fprintf(qfile, "  -sf '%s' \\\n", a_header->remotename);
  fprintf(qfile, "    >& %s\n\n", resfilename);

  /* mail the MEME results */
  fprintf(qfile, "%s/mailer\\\n", meme_bin);
  fprintf(qfile, "  %s\\\n", a_header->address);
  fprintf(qfile, "  \'MEME job %ld results: %s\' %s %s\n\n",
    job, a_header->description, resfilename, type);

  /*
   run MAST on the MEME output and mail results back
  */
  fprintf(qfile, "$bin/make_logodds %s /dev/null > /dev/null\n", 
    resfilename);
  fprintf(qfile, "if ($status == 0) then\n");
  fprintf(qfile, "  set mast_out = $dir/mast.%ld.results\n", job);
  fprintf(qfile, "  $bin/mast %s %s -mf \'<motifs found in %s>\'\\\n",
    resfilename, meme_input, a_header->remotename);
  fprintf(qfile, "    -df \'%s\' -stdout %s > $mast_out\n",
    a_header->remotename, text);
  fprintf(qfile, "  %s/mailer\\\n", meme_bin);
  fprintf(qfile, "    %s\\\n", a_header->address);
  fprintf(qfile, "    \'MAST job %ld results: %s\' $mast_out %s\n",
    job, a_header->description, type);
  fprintf(qfile, "  /bin/rm $mast_out\n");
  fprintf(qfile, "endif\n\n");

  /* save the final time */
  fprintf(qfile, "set t2 = `date %s`\n", dfmt);

  /* log this job */
  fprintf(qfile, "touch $dir/meme-log\n");		/* make sure log-file exists */
  fprintf(qfile, 
    "echo `hostname` %9ld submit: `cat %s` start: $t1 end: $t2 %s %s %s >> $dir/meme-log\n\n", 
    job, datefilename, a_header->switches, a_header->loginfo, a_header->address
  );

  /* fprintf(qfile, "exit 0 \n" );	*/

  /* delete the files we have created */
  fprintf(qfile, "/bin/rm -f %s\n", meme_input);	/* sequence file */
  fprintf(qfile, "/bin/rm -f %s\n", datefilename);	/* submission date file */
  fprintf(qfile, "/bin/rm -f %s\n", resfilename);	/* output of MEME */
  fprintf(qfile, "/bin/rm -f %s\n", qfilename);	/* this script */

  /* Close the script file. */
  fclose(qfile);

  /* free the file names created here */
  myfree(resfilename);
  myfree(datefilename);
}

/***********************************************************************
 *
 * void confirm
 *
 * Send a confirmation message to the user
 *
 ***********************************************************************/
void confirm( HEADER *a_header, char * meme_bin )
{
  char command[2024];                   /* confirmation message command */
  int job = getpid();                   /* get the job number */
  int result;
 
  /* make the message */
  sprintf(command, 
"%s/mailer %s \'MEME job %d confirmation: \' -stdin << END\n"
"Your MEME search request %d is being processed.\n"
"You should receive two subsequent messages containing:\n"
"	1) MEME search results\n"
"	2) MAST analysis of your sequences using the MEME results\n"
"END\n",
    meme_bin, a_header->address, job, job);
 
  /* send it */
  if ((result = system(command)) != SYSTEM_OK) {
    fprintf(stderr, "Error level %d from \'%s\'.\n", result, command);
    fflush(stderr);
    exit(30);
  }
}
 
/***********************************************************************
 *
 * void submit_job
 *
 * Submit the given job file to queueing system.
 *
 ***********************************************************************/
void submit_job(char *qfilename)
{
  int result;
  char command[2000]; 

  sprintf(command, QCMD, qfilename);
  if ((result = system(command)) != SYSTEM_OK) {
    fprintf(stderr, "Error level %d from \'%s\'.\n", result, command);
    fflush(stderr);
    exit(30);
  }
}


/***********************************************************************
 *
 * doit
 *
 * Do whatever the server is supposed to do with incoming data.
 *
 ***********************************************************************/
void doit (int sock)
{
  char *meme_input;
  char *qfilename;
  HEADER a_header;
  static char *meme_logs = NULL;     /* meme logs directory */
  static char *meme_bin = NULL;     /* meme bin directory */

  meme_logs = getenv("MEME_LOGS");
  if (meme_logs == NULL) {
    fprintf(stderr, "You must define environment variable MEME_LOGS\n");
    exit(1);
  } 

  meme_bin = getenv("MEME_BIN");
  if (meme_bin == NULL) {
    fprintf(stderr, "You must define environment variable MEME_BIN\n");
    exit(1);
  }


  /* get unique filenames */
  meme_input = newfilename("data", meme_logs);
  qfilename = newfilename("script", meme_logs);

  /* Get the file from the client. */
  receive_file(sock, meme_input, &a_header);
  send_ack(sock);

  /* Write out an NQS job file. */
  make_q_script(&a_header, meme_input, qfilename, meme_bin, meme_logs);
  
  /* Submit the job file to the queue. */
  submit_job(qfilename);

  /* Send confirmation message */
  confirm(&a_header, meme_bin);
 
  /* free filenames */
  myfree(meme_input);
  myfree(qfilename);
}

/***********************************************************************
  void reapchild


  SIGCHLD handler. Supplied by Mike Wan at SDSC. (11/9/95)

***********************************************************************/
#include <sys/wait.h>
#include <signal.h>
static void reapchild (int sig)
{
  int status;
  pid_t waitval;
 
  while ( (waitval=waitpid (-1, &status, WNOHANG)) > 0) {
    if (sig == SIGCHLD) {
       /*printf("Child process %d exited.\n", (int) waitval);*/
    } else {
      printf( "Received signal %d from process %d\n", sig, (int) waitval );
    }
    fflush(stdout);
    signal (SIGCHLD, reapchild);
  }
}

/***********************************************************************
 *
 * main
 *
 ***********************************************************************/
extern int main(
  int argc,
  char *argv[]
)
{
  #if defined(ibmrs6000)
  size_t len;
  #else
  int len;
  #endif
  int sock, new_socket;
  struct sockaddr_in acceptor;
  struct sockaddr_in connector;

  if (argc < 4) {
    /*
    printf("Usage: meme-server <socket> <meme> <qcmd>\n");
    printf("	<socket>  		socket number to listen on\n");
    */
    printf("	<command>		command for server to execute\n");
    printf("	<qcmd>			queueing command including %%s\n");
    fflush(stdout);
    exit(1);
  }

  /* get root directory for meme system */
  SOCKET_NUMBER = atoi(argv[1]);
  MEME = argv[2];
  QCMD = argv[3];

  printf("MEME server initialized.\n");
  printf("command is: %s\n", MEME);
  printf("qcmd is: %s\n", QCMD);
  fflush(stdout);

  /* Parent needs to handle the SIGCHLD signal before the
     child can exit properly. */
  signal (SIGCHLD, reapchild);

  len = sizeof(acceptor);
  memset(&acceptor, 0, len);
  acceptor.sin_family = AF_INET;
  acceptor.sin_port = 0xffff & (htons(SOCKET_NUMBER));
  
  connector.sin_addr.s_addr = INADDR_ANY;

  /* Create a new socket. */
  sock = socket(AF_INET, SOCK_STREAM, 0);

  /* Bind the socket. */
  if (bind(sock, (struct sockaddr *)&acceptor, len) < 0) {
    perror("Node bind error");
    fprintf(stderr, "Node errno is %d\n", errno);
    fflush(stderr);
    exit(14);
  }

  /* Listen on the socket. */
  if (listen(sock, 5) < 0 ) {
    perror("Node listen");
    fflush(stderr);
    exit(15);
  }
  
  /* Loop forever, forking for each incoming request. */
  for(;;)
    {
      /* Block until a request comes in, then create a new socket. */
      new_socket = accept(sock, (struct sockaddr *)&acceptor, &len);
      if (new_socket < 0) {
        /* allow interrupt caused by child process terminating */
        if (errno == EINTR) {
           continue;            /* go to top of loop */
        } else {
          perror("Node accept");
          fflush(stderr);
          exit(16);
        }
      }
 
      /* On a fork, the parent gets the child's process ID back. */
      if (fork()==0) {     /* This means I am child. */
	close(sock);       /* Ignore the old socket. */
	doit(new_socket);  /* <- Primary server activity here. */
	exit(0);           /* After execution, the child exits. */
      }

      /* The parent closes the new socket. */
      close(new_socket);
    }

  /* NOTREACHED */
  return(0);
}
