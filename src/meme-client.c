/*
 * $Id: meme-client.c,v 1.2 2006/01/03 06:37:00 tbailey Exp $
 * 
 * $Log: meme-client.c,v $
 * Revision 1.2  2006/01/03 06:37:00  tbailey
 * Fix "MEME_BINbin" bug in mast-client.txt.
 * Fix indentation in meme-server.c and meme-client.c.
 *
 * Revision 1.1.1.1  2005/07/29 17:21:42  nadya
 * Importing from meme-3.0.14, and adding configure/make
 *
 */

/*
 AUTHOR: Bill Grundy
 CREATE DATE: 07-06-95
 USAGE: meme-client <socket> <hostname> <filename>

 DESCRIPTION:

   This program is part of the MEME WWW interface. It is one half of
the client-server interface between the http daemon (the client) and
the Paragon server.

   'meme-client' gets called whenever a user submits data. The program
opens a socket to the server and sends the parameters and data across
the socket.

   'meme-client' does not wait for the Paragon server to process the
data. As soon as the data is delivered, the client process exits.  The
Paragon server e-mails the results to the user.

  The program is based upon code supplied by Allan Snavely at SDSC.


Error messages:
---------------
  10 : Wrong number of command-line arguments
  11 : Host name not found
  12 : Cannot connect
  13 : Error opening file
  14 : Error writing to socket

*/

/***********************************************************************
 * Include files
 ***********************************************************************/
#include "meme-cs.h"

/***********************************************************************
 *
 * void get_cmdline
 *
 ***********************************************************************/
void get_cmdline(
  int argc, 
  char *argv[], 
  int *socket,			/* socket number */
  char *servername, 		/* name of cpu running server */
  char *datafilename,		/* name of data file */
  int *ntries 			/* number of times to try connection */
)
{

  /* Make sure we've got the right number of arguments. */
  if (argc < 4 || argc > 5) {
    /*
    fprintf(stderr, 
      "USAGE: meme-client <socket> <hostname> <data filename> [<ntries>]\n");
    */
    exit(10);
  }

  /* Grab the args. */
  *socket = atoi(argv[1]);
  strcpy(servername, argv[2]);
  strcpy(datafilename, argv[3]);
  if (argc == 5) {
    *ntries = atoi(argv[4]);
  }
}

/***********************************************************************
 * 
 * int daemen_connect
 *
 * Initiate a connection between the current machine and a specified
 * server.
 *
 ***********************************************************************/
int daemen_connect(char *serv, int socknum, int ntries)
{
  int sock = 0, len, tries;
  struct sockaddr_in connector;
  struct hostent *hp;

  len = sizeof(connector);
  memset(&connector,0,len);	/* clear connector struct */
  connector.sin_family = AF_INET;
  connector.sin_port   = htons(socknum);

  /* Get the name of the server. */
  hp = gethostbyname(serv);

  /* Indicate error if host name not found. */
  if(!hp) {
    fprintf(stderr, "Host %s not found\n", serv);
    exit(11);
  }

  memcpy(&connector.sin_addr, hp->h_addr, hp->h_length);

  /* Attempt a connection. */
  for (tries = 0; tries < ntries; tries++) {
    sock = socket(AF_INET, SOCK_STREAM, 0);
    if (connect(sock,(struct sockaddr *)&connector,len) >= 0) {
      break;				/* suceeded */
    } else {				
      if (errno != ECONNREFUSED) {	/* failed */
	perror("host connect error");
      } else {
	printf("host connect failed\n");
      }
      if (tries == ntries - 1) {
        fprintf(stderr, "Giving up trying to connect to %s\n", serv);
        exit(12);
      }
      fprintf(stderr, "retrying...\n");
      close(sock);
      sleep(1);
    } /* ntries */
  }
  return(sock);

}

/***********************************************************************
 * 
 * meme_send_file
 *
 * Send a file across a socket.
 *
 ***********************************************************************/
void meme_send_file(int sock, char *filename)
{
  int filedescriptor;
  char achar;

  /* Open the file in read-only mode. */
  if ((filedescriptor = open(filename, 0)) == -1) {
    fprintf(stderr, "Error opening file %s\n", filename);
    close(sock);
    exit(13);
  }
  
  /* Transfer the file byte by byte. */
  while (read(filedescriptor, &achar, 1) == 1) {
    /* fprintf(stderr, "%c\n", achar);*/
    if (writen(sock, &achar, 1) != 1) {
      fprintf(stderr, "Error writing to socket.\n");
      close(sock);
      exit(14);
    }
  }

  /* Close the socket and the file. */
  close(sock);
  close(filedescriptor);
}


/***********************************************************************
 * 
 * main
 *
 ***********************************************************************/
int main(int argc, char *argv[])
{
  char servername[80];
  char datafilename[80];
  int socket;			/* integer format socket number */
  int sock;			/* htons format socket number */
  int ntries = 50;		/* number of times to try connection */

  /* Get the command line arguments. */
  get_cmdline(argc, argv, &socket, servername, datafilename, &ntries);

  /* Connect to the named server. */
  sock = daemen_connect(servername, socket, ntries);
  /*fprintf(stdout, "Connected to %s\n", servername);*/

  /* Send the file across the socket. */
  meme_send_file(sock, datafilename);

  return(0);
}
