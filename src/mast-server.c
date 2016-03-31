/*
 * $Id: mast-server.c,v 1.6 2005/10/28 22:44:19 nadya Exp $
 * 
 * $Log: mast-server.c,v $
 * Revision 1.6  2005/10/28 22:44:19  nadya
 * fix typo in "rm"
 *
 * Revision 1.5  2005/09/15 21:27:12  nadya
 * fix db path for uploaded files
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
 * Revision 1.1.1.1  2005/07/29 17:16:51  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*
	mast-server

	Usage: mast-server <socket> <command>
		<socket>  socket number to listen on
		<command> command to execute; quote if spaces

        If the file sent to the socket contains just the word "ping",
        the program exits with success.
 
	6/02/96; created; tbailey

Format for input file read from the socket:
----------------------
ADDRESS <return address> 
PROGRAM mast
DESCRIPTION <1 line of text>
MOTIFS <original name of motif file>
DB_TYPE local | uploaded
DATABASE <name of database file to search (../mast_databases/ will be prepended)>
ALPHABET <alphabet used in database>
SWITCHES <other switches>
LOGINFO <information to log>
BEGIN_MOTIFS
<logodds matrices>
[
%
BEGIN_DB
<uploaded FASTA sequences>
%
]

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
  22 : Illegal number of motifs
*/

/***********************************************************************
 * Include files
 ***********************************************************************/
#include "meme-cs.h"
#include <signal.h>
#include <sys/wait.h>
#include <sys/errno.h>

/***********************************************************************
 * Constants
 ***********************************************************************/
#define MAXHEADERLINE 256
#define NPROCS 1

/*
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
  22 : Illegal number of motifs
*/

/***********************************************************************
 * Include files
 ***********************************************************************/
#include "meme-cs.h"
#include <signal.h>
#include <sys/wait.h>
#include <sys/errno.h>

/***********************************************************************
 * Constants
 ***********************************************************************/
#define MAXHEADERLINE 256
#define NPROCS 1
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
  char description[MAXHEADERLINE];
  char motif_file[MAXHEADERLINE];
  char db_type[MAXHEADERLINE];
  char database[MAXHEADERLINE];
  char alphabet[MAXHEADERLINE];
  char switches[MAXHEADERLINE];
  char loginfo[MAXHEADERLINE];
} HEADER;


/***********************************************************************
 * Globals 
 ***********************************************************************/

int SOCKET_NUMBER;		/* socket to listen on */
char *MAST;			/* name of mast executable */


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
    a_string[index] = tolower((int) a_string[index]);
}

/***********************************************************************
 *
 * void read_header_item
 *
 * Read one header line from a socket into a HEADER data structure.
 *
 ***********************************************************************/
void read_header_item(HEADER *a_header, int sock, char *fieldname)
{
  int i;
  char buffer[MAXHEADERLINE];
  char *info;
  
  /* Read the line into a local buffer. */
  if (readline(sock, buffer, MAXHEADERLINE) == 0) {
    fprintf(stderr, "Error reading from socket.\n");
    fflush(stderr);
    exit(11);
  }

  /* Check to see if this is a ping */
  if (strcmp(buffer, "ping") == 0) {
    exit(0);                            /* ping; exit with OK status */
  }
 
  /* Search for the field name in the line. */
  if ((info = (char *)strstr(buffer, fieldname)) == NULL) {
    fprintf(stderr, "Field %s not found in %s.\n", fieldname, buffer);
    fflush(stderr);
    exit(20);
  }

  /* Move past the end of the field name. */
  info += strlen(fieldname);
  if (*info != '\0') info++;

  /* Get the first character of the field name. */
  switch (fieldname[0]) {
    case 'A' : 
      switch (fieldname[1]) {
	case 'D' :  /* ADDRESS */
	  strcpy(a_header->address, info);
	  break;
	case 'L' :  /* ALPHABET */
	  strcpy(a_header->alphabet, info);
	  break;
      }
      break;
    case 'M' :  /* MOTIFS */
      strcpy(a_header->motif_file, info);
      break;
    case 'P' :  /* PROGRAM */
      strcpy(a_header->program, info);
      break;
    case 'D' : 
      switch (fieldname[1]) {
	case 'E' :  /* DESCRIPTION */
	  strcpy(a_header->description, info);
	  break;
	case 'A' :  /* DATABASE */
	  strcpy(a_header->database, info);
	  break;
	case 'B' :  /* DB_TYPE */
	  strcpy(a_header->db_type, info);
	  break;
      }
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
    case 'B' :  /* BEGIN_LO or BEGIN_MOTIFS */
      /* Don't do anything except skip header line */
      break;
    default: 
      fprintf(stderr, "Invalid fieldname: %s\n", fieldname);
      fflush(stderr);
      exit(21);
  }
} /* read_header_item */


/***********************************************************************
 *
 * write_data_file
 *
 * Read in some lines from a socket and write them to a file.
 *
 ***********************************************************************/
void write_data_file(int sock, char *filename, char terminator)
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
      
  /* Read from the socket until line containing only terminator reached. */
  for (nread = readn(sock, &achar, 1); nread == 1 && achar != terminator;
       nread = readn(sock, &achar, 1)) {
    
    /*printf("%c", achar);*/

    /* Echo incoming data to the file. */
    if (write(filedescriptor, &achar, 1) != 1) {
      fprintf(stderr, "Error writing to file %s.\n", filename);
      fflush(stderr);
      exit(12);
    }
  }

  /* Throw away newline after terminator. */
  if (achar == terminator) {
    nread = readn(sock, &achar, 1);
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
void receive_file(
  int sock, 
  char *motif_filename, 
  char *db_filename,
  HEADER *a_header
)
{
  /* Read in the header. */
  read_header_item(a_header, sock, "ADDRESS");
  read_header_item(a_header, sock, "PROGRAM");
  read_header_item(a_header, sock, "DESCRIPTION");
  read_header_item(a_header, sock, "MOTIFS");
  read_header_item(a_header, sock, "DB_TYPE");
  read_header_item(a_header, sock, "DATABASE");
  read_header_item(a_header, sock, "ALPHABET");
  read_header_item(a_header, sock, "SWITCHES");
  read_header_item(a_header, sock, "LOGINFO");
  read_header_item(a_header, sock, "BEGIN_MOTIFS");
  /* Echo the motifs to a file. */
  write_data_file(sock, motif_filename, '%');
  if (strcmp(a_header->db_type, "uploaded") == 0) {
    /* Echo the uploaded sequence data to a file. */
    read_header_item(a_header, sock, "BEGIN_DB");
    write_data_file(sock, db_filename, '%');
  }
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
 *	mast.<pid>.<filename>
 *
 ***********************************************************************/
char *newfilename(char *filename, char * dir)
{
  int length = 6 + 20 + strlen(filename) + 1 + strlen(dir);
  char *name = (char *) mymalloc(length);

  sprintf(name, "%s/mast.%ld.%s", dir, (long) getpid(), filename);
  return name;
}

/***********************************************************************
 *
 * void make_q_script
 *
 * Create a script to run the job and mail back the results.
 *
 ***********************************************************************/
void make_q_script(HEADER *a_header, char *motiffilename, 
                   char *dbfilename, char *qfilename, char *rfilename, 
                   char *meme_bin, char *meme_logs, char *meme_db)
{
  FILE *qfile;
  int job = getpid();			/* get the job number */
  char *type;				/* type of file to mail */

  /* Open the script file */
  if ((qfile = fopen(qfilename, "w")) == NULL) {
    fprintf(stderr, "Error opening file %s.\n", qfilename);
    fflush(stderr);
    exit(10);
  }

  /* put in script header */
  fprintf(qfile, "#!/bin/csh \n\n");
  fprintf(qfile, "# This shell script was automatically generated ");
  fprintf(qfile, "by the MAST server (mast-server.c).\n\n");
  fprintf(qfile, "set dir = %s\n", meme_logs);
  fprintf(qfile, "set bin = %s\n\n", meme_bin);

  /* log this job */
  fprintf(qfile, "touch $dir/mast-log\n");           /* make sure log-file exists */
  fprintf(qfile, "set db = '%s'\n", a_header->database);
  fprintf(qfile, 
    "echo `hostname` %9d `date` \'%s\' $db %s %s %s %s >> $dir/mast-log\n\n",
    job, a_header->motif_file, a_header->alphabet, 
    a_header->switches, a_header->loginfo, a_header->address 
  );

  /* run mast and mail the results back */
  type = strstr(a_header->description, "(Use web browser to view results)") ? 
    "-html" : "";
  if (strcmp(a_header->db_type, "uploaded") == 0) {
    fprintf(qfile, "set db = '%s'\n", dbfilename);
  } else {
    fprintf(qfile, "set db = %s/%s\n", meme_db, a_header->database);
  }
  fprintf(qfile, "set bfile = ''\n");
  if (strstr(a_header->switches, "-dna")) {	/* translating dna */
    fprintf(qfile, "if (-e $db.xbfile) set bfile = \"-bfile $db.xbfile\"\n");
  } else {					/* not tranlating dna */
    fprintf(qfile, "if (-e $db.bfile) set bfile = \"-bfile $db.bfile\"\n");
  }
  fprintf(qfile, "set nseqs = ''\n");
  fprintf(qfile, "if (-e $db.nseqs) set nseqs = \"-minseqs `cat $db.nseqs`\" \n");
  fprintf(qfile, "cat $db | %s \\\n", MAST);
  fprintf(qfile, "  %s -mf '%s' -stdin -remcorr \\\n", motiffilename, a_header->motif_file);
  fprintf(qfile, "  %s \\\n", a_header->switches);	/* overridable */
  fprintf(qfile, "  $nseqs $bfile -a %s \\\n", a_header->alphabet);
  fprintf(qfile, "  -stdout -nostatus >& %s \n", rfilename);
  fprintf(qfile, "%s/mailer\\\n",  meme_bin);
  fprintf(qfile, "  %s\\\n", a_header->address);
  fprintf(qfile, "  \'MAST job %d results: %s\' %s %s\n\n", 
    job, a_header->description, rfilename, type);
  /* fprintf(qfile, "exit 0 \n"); */

  fprintf(qfile, "/bin/rm -f '%s' \n", motiffilename);
  fprintf(qfile, "/bin/rm -f '%s' \n", dbfilename);
  fprintf(qfile, "/bin/rm -f '%s' \n", rfilename);
  fprintf(qfile, "/bin/rm -f '%s' \n", qfilename);


  /* Close the job file. */
  fclose(qfile);
}

/***********************************************************************
 *
 * void confirm
 *
 * Send a confirmation message to the user
 *
 ***********************************************************************/
void confirm( HEADER *a_header, char * meme_bin)
{
  char command[2024];			/* confirmation message */
  int job = getpid();			/* get the job number */
  int result;

  /* make the message */
  sprintf(command, 
"%s/mailer %s \'MAST job %d confirmation: \' -stdin << END\n"
"Your MAST search request %d is being processed:\n"
"  Motif file: %s\n"
"  Database to search: %s\n"
"You should receive a subsequent message containing:\n"
"       1) MAST search results\n"

"END\n",
    meme_bin, a_header->address, job, 
    job, a_header->motif_file, a_header->database);

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
 * Submit the given job file to NQS or the shell.
 *
 ***********************************************************************/
void submit_job(char *qfilename)
{
  char command[1024]; 
  int result;

  /* set up the batch command depending on the type of OS */
  sprintf(command, "/bin/csh %s &", qfilename);

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
  char *motiffilename;
  char *dbfilename;
  char *qfilename;
  char *rfilename;
  HEADER a_header;
  static char *meme_logs = NULL;    /* meme logs directory */
  static char *meme_bin = NULL;     /* meme bin directory */
  static char *meme_db = NULL;      /* meme databases directory */

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

  meme_db = getenv("MEME_DB");
  if (meme_db == NULL) {
      fprintf(stderr, "You must define environment variable MEME_DB\n");
      exit(1);
  }


  /* get unique filenames */
  motiffilename = newfilename("motifs", meme_logs);
  dbfilename = newfilename("uploaded_db", meme_logs);
  qfilename = newfilename("script", meme_logs);
  rfilename = newfilename("results", meme_logs);

  /* Get the file from the client. */
  receive_file(sock, motiffilename, dbfilename, &a_header);
  send_ack(sock);

  /* Write out an NQS job file. */
  make_q_script(&a_header, motiffilename, dbfilename, qfilename, rfilename,
                meme_bin, meme_logs, meme_db);

  /* Log the job. */
  printf("\njob: %d\n", (int) getpid());
  fflush(stdout);
  
  /* Submit the job file to the queue. */
  submit_job(qfilename);

  /* Send confirmation message */
  confirm(&a_header, meme_bin);

  /* free filenames */
  myfree(motiffilename);
  myfree(qfilename);
  myfree(dbfilename);
  myfree(rfilename);
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
  int len, sock, new_socket;
  struct sockaddr_in acceptor;
  struct sockaddr_in connector;


  if (argc < 3) {
    /*
    fprintf(stderr, "Usage: mast-server <socket> <command>\n");
    fprintf(stderr, "	<socket>  	socket number to listen on\n");
    */
    fprintf(stderr, "	<command> 	command to for server to execute\n");
    exit(1);
  }

  /* get root directory for meme system */
  SOCKET_NUMBER = atoi(argv[1]);
  MAST = argv[2];

  system("echo MAST server initialized on CPU: `hostname`");
  printf("command is: %s\n", MAST);
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
           continue;		/* go to top of loop */
        } else {
	  perror("Node accept");
          fflush(stderr);
	  exit(16);
        }
      }

      /* Create a child process to perform the command */
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
  exit(0);
}
