## $Id: configure.ac,v 1.30.4.5 2006/01/31 20:53:31 nadya Exp $
##
## Process this file with autoconf to produce a configure script.
##
## # $Log: configure.ac,v $
## # Revision 1.30.4.5  2006/01/31 20:53:31  nadya
## # use spaces to delimit file names (for cygwin)
## #
## # Revision 1.30.4.4  2006/01/31 20:02:52  nadya
## # use ":" delimiter in file list for cygwin
## #
## # Revision 1.30.4.3  2006/01/31 19:06:05  nadya
## # change CYGWIN to Cygwin. The former has a special meaning to OS.
## # add backslashes to the file names to prevent newline translation
## # problems on cygwin.
## #
## # Revision 1.30.4.2  2006/01/24 02:24:29  nadya
## # EXT is obsoleted by using printf()
## #
## # Revision 1.30.4.1  2006/01/24 02:16:01  nadya
## # add EXT for echo output (by platform)
## #
## # Revision 1.30  2005/12/16 23:31:55  nadya
## # if parallel is disabled make server disabled.
## # currently this is for mac and more likely cygwin platforms
## # or when no mpi is installed.
## #
## # Revision 1.29  2005/12/16 05:08:35  tbailey
## # Add support for CYGWIN operating system.
## #
## # Revision 1.28  2005/12/13 01:01:37  nadya
## # add WANT_STARTUP for making platform dependent scripts in etc/
## #
## # Revision 1.27  2005/10/28 23:53:21  nadya
## # make location of LOGS configurable
## #
## # Revision 1.26  2005/10/25 21:25:20  nadya
## # make LOGS/ location configurable during "configure" step
## #
## # Revision 1.25  2005/10/25 18:51:12  nadya
## # rm debug output
## #
## # Revision 1.24  2005/10/20 00:20:51  tbailey
## # *** empty log message ***
## #
## # Revision 1.23  2005/10/19 19:36:23  nadya
## # for irix v. 6*: redefine host and rm lmpi from LIBS.
## #
## # Revision 1.22  2005/10/14 21:20:49  nadya
## # add --with-mpicc switch for providing a specific mpicc.
## # add mpcc to MPICC for automatic detection
## #
## # Revision 1.21  2005/10/07 19:23:27  nadya
## # check uid of the user who runs configure.
## # set MEME_USER to uid if --with-user is not provided.
## #
## # Revision 1.20  2005/10/06 22:04:50  nadya
## # add -O flags when compiling on sun
## # add -lm when linking on sun
## # upgrade distro version
## #
## # Revision 1.19  2005/10/04 23:47:50  nadya
## # lower prerequisite autoconf verison to 2.53
## #
## # Revision 1.18  2005/10/02 06:11:27  nadya
## # test for qstat, if not found disable parallel compilation
## # add switch for compiling in serial mode only
## # set no web installation as a default
## # clean all options comments
## #
## # Revision 1.17  2005/09/21 20:38:44  nadya
## # add user to run servers
## #
## # Revision 1.16  2005/09/15 21:28:31  nadya
## # Wilfred wanted 0 in verison
## #
## # Revision 1.15  2005/09/13 22:07:24  nadya
## # move gnu-related files in config/
## #
## # Revision 1.14  2005/09/13 19:04:24  nadya
## # set CONTACT to meme@host, othwerwise email from the compute nodes has wrong email address.
## # update help strings
## # set version to 3.5
## #
## # Revision 1.13  2005/09/01 18:47:32  nadya
## # add missing slash in url
## #
## # Revision 1.12  2005/09/01 01:31:17  nadya
## # change check to test, otherwise all dirs are traversed automatically.
## #
## # Revision 1.11  2005/09/01 00:52:09  nadya
## # put "no" in rught places for server and client
## #
## # Revision 1.10  2005/08/31 22:49:58  nadya
## # really do these changes that were claimed in v.1,8
## #
## # Revision 1.9  2005/08/31 22:10:57  nadya
## # update configuration output
## # add info on running make
## #
## # Revision 1.8  2005/08/31 21:43:46  nadya
## # set default prefix to $HOME
## # add check for ld
## # fix list of math functions
## # fix prefix in MPI_CMD, MEME_DIRECTORY
## # add LAMBIN for use in full path in scripts
## # add enable-client option
## # enable server and client build by defaul
## # remove wc from lampi test
## #
## # Revision 1.7  2005/08/25 03:20:34  nadya
## # change test for SGE: remove "wc" from the test.
## # somw versions of wc pad extra spaces to the output.
## #
## # Revision 1.6  2005/08/24 19:27:53  nadya
## # add MACH_ID for using ini meme.csh
## # add MEME_DB for database location
## #
## # Revision 1.5  2005/08/12 18:01:23  nadya
## # add variables for running meme
## #
## # Revision 1.4  2005/08/09 23:08:49  nadya
## # add sed
## #
## # Revision 1.3  2005/08/01 01:56:21  nadya
## # rm msg line
## #
## # Revision 1.2  2005/08/01 01:36:26  nadya
## # change location of WEBDIR to depend on prefix
## #
## # Revision 1.1.1.1  2005/07/30 02:19:40  nadya
## # Importing from meme-3.0.14, and adding configure/make
## #
##

AC_PREREQ(2.53)
AC_INIT(src/alphtype.c)
AM_CONFIG_HEADER([config.h])
AC_CONFIG_AUX_DIR(config)
AM_INIT_AUTOMAKE(memeris, 1.0)
AC_PREFIX_DEFAULT(${HOME})

# Checks for programs.
AC_PROG_CC
AC_PROG_RANLIB
AC_CANONICAL_HOST
AC_PROG_INSTALL
AC_PROG_MAKE_SET
AC_PATH_PROG(LD, ld)
AC_PATH_PROG(SED, sed)
AC_PATH_PROG(MPICC, mpicc mpcc)
AC_PATH_PROG(MAIL, mail)
AC_PATH_PROG(SENDMAIL, sendmail, [/usr/lib/sendmail], [$PATH:/usr/lib:/usr/sbin])
AC_PATH_PROG(PERL, perl)
AC_PATH_PROG(QSTAT, qstat)

AC_DEFINE_UNQUOTED(HOSTNAME, "`hostname`", [Build host name])

# Checks for libraries.
AC_CHECK_LIB([m], 
             [isnan, sqrt, log, exp, fabs, pow, floor, ceil], 
             [], 
             [AC_MSG_WARN([Math library not found])])

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_SYS_WAIT
AC_CHECK_HEADERS([arpa/inet.h fcntl.h limits.h netdb.h netinet/in.h stdlib.h string.h strings.h sys/socket.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_HEADER_STDBOOL
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_TYPE_SIGNAL
AC_FUNC_STAT
AC_CHECK_FUNCS([floor getcwd gethostbyname isascii memset pow rint socket sqrt strchr strcspn strdup strspn strstr])
AC_CHECK_FUNCS([strlcpy], [have_strlcpy=yes], [have_strlcpy=no])
AM_CONDITIONAL(WANT_STRLCPY, test ["$have_strlcpy" = no])

# Check for compiler-specific settings
ac_mpicc_flags="-O"
if test "${CC}" = "gcc"; then
  ac_cc_opt="-O3"
  ac_cc_warn="-Wall"
  CFLAGS="$ac_cc_warn -DUNIX -D__USE_FIXED_PROTOTYPES__"
else
  ac_cc_opt="-O"
  ac_cc_warn=""
  CFLAGS="$ac_cc_warn -DUNIX -D__USE_FIXED_PROTOTYPES__"
fi

# who is installing 
USER=`whoami`

# Check for platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
	  have_startup="yes"
      ;;
   *aix*)
      ac_mpicc_flags="-O3 -qstrict -qarch=pwr3 -qtune=pwr3"
      AC_DEFINE(AIX, [], [Using AIX])
      MACHID="aix"
      if /usr/sbin/lsattr -EHl proc0 | grep POWER >/dev/null 2>&1; then
          MACHID=rs6000
      fi
      ;;
   *cray*)
      ac_mpicc_flags="-O3 -h msglevel_3"
      AC_DEFINE(cray3e, [], [Using cray3e])
      MACHID="cray"
      ;;
   *irix6*)
      ac_mpicc_flags="-O3"
      AC_DEFINE(IRIX64, [], [Using IRIX64])
      MACHID="irix"
	  have_startup="yes"
      ;;
   *irix*)
      ac_mpicc_flags="-O3"
      MACHID="irix"
	  have_startup="yes"
      ;;
   *solaris*|*sunos*)
      if test "${CC}" != "gcc"; then
          ac_cc_opt="-O -Xa"
      fi
      AC_DEFINE(SunOS, [], [Using SunOS])
      AC_CHECK_LIB([nsl], [gethostbyname], [], [AC_MSG_ERROR([Nsl library not found])])
      AC_CHECK_LIB([socket], [socket], [], [AC_MSG_ERROR([Socket library not found])])
      LIBS="-lsocket -lnsl -lm"
      MACHID="solaris"
	  have_startup="yes"
      ;;
esac

AM_CONDITIONAL(WANT_STARTUP, test ["$have_startup" = yes])

####################################
# Check for configuration switches
####################################

# debugging
AC_ARG_ENABLE(
    debug,
    [  --enable-debug          Produce an executable with debugging symbols. Default is disabled.],
    [if test "$enableval" = "yes"; then
        DEBUG_FLAGS="-g";
     fi],
    [DEBUG_FLAGS=""])

# optimization
AC_ARG_ENABLE(
    opt,
    [  --enable-opt            Produce an executable with optimization. Default compiler dependent.],
    [if test "$enableval" = "yes"; then
       OPT_FLAGS="$ac_cc_opt";
     fi],
    [OPT_FLAGS="$ac_cc_opt"])

# queue type
AC_ARG_WITH(
    queue,
    [  --with-queue=QUEUE      QUEUE type (SGE or PBS). Default is the type of found qstat.],
    [QTYPE="$withval"],
    [QTYPE=""])

if test "$QTYPE" = "" -a "$QSTAT" != "" ; then
   info=`qstat -help | grep SGE`;
   if test -n "$info" ; then
      QTYPE="sge";
   else
      QTYPE="pbs";
   fi
fi

# get mpi info, and check if parallel version can be compiled
AC_ARG_WITH(
    mpicc,
    [  --with-mpicc=MPICC      Set MPICC to use. Default is found mpicc.],
    [MPICC="$withval"],
    [MPICC="$MPICC"])
    
AC_ARG_WITH(
    mpidir,
    [  --with-mpidir=MPIDIR    Set MPI installation path. Default is path of found mpicc.],
    [mpidir="$withval";
     MPICC=$mpidir/bin/mpicc;
     MPILIBDIR=$mpidir/lib;
     MPIINC="-DMPI -I$mpidir/include";
     AC_CHECK_LIB([mpi], [MPI_Comm_rank])
     parallel="yes";],
    [if test "$MPICC" = ""; then
        MPILIBDIR="";
        MPIINC="";
        parallel="no";
     else
        mpidir=`echo $MPICC | sed -e 's/\/bin\/.*//'`
        MPILIBDIR="$mpidir/lib"
        MPIINC="-DMPI -I$mpidir/include"
        AC_CHECK_LIB([mpi], [MPI_Comm_rank])
        parallel="yes";
     fi])

if test "$QSTAT" = ""; then
        parallel="no";
fi

# enable only serial installation
AC_ARG_ENABLE(
    serial,
    [  --enable-serial         Produce only serial executables. Default parallel and serial.],
    serial=$enableval, serial=no)
if test "$serial" = "yes"; then
    parallel="no";
fi

if test "$parallel" = "yes"; then
   lammpi=`strings $MPICC | grep liblam`;
   if test -n "$lammpi" ; then
      MPI_CMD="\${prefix}/bin/lrun -n";
      MPITYPE="lam"
      LAMBIN="$mpidir/bin"
   else
      MPI_CMD="$mpidir/bin/mpirun -np";
      MPITYPE="mpich"
   fi
fi
AM_CONDITIONAL(PARALLEL, test ["$parallel" != no])

# enable server installation
AC_ARG_ENABLE(
    server,
    [  --enable-server         Produce server executables. Default is enabled.],
    server=$enableval, server=yes)

if test "$parallel" = "no"; then
    server="no";
fi
AM_CONDITIONAL(SERVER, test ["$server" != no])

# enable client installation
AC_ARG_ENABLE(
    client,
    [  --enable-client         Produce client executables. Default is enabled.],
    client=$enableval, client=yes)
AM_CONDITIONAL(CLIENT, test ["$client" != no])

# enable website installation
AC_ARG_ENABLE(
    web,
    [  --enable-web            Build website applications. Default is disabled.],
    web=$enableval, web=no)
AM_CONDITIONAL(WEBSITE, test ["$web" != no])

# set website url
AC_ARG_WITH(
    url,
    [  --with-url=URL          Set url for the website. Default http://`hostname`/meme.],
    [URL=$withval],
    [URL="http://$HOSTNAME/meme"])

# set meme user 
AC_ARG_WITH(
    user,
    [  --with-user=USER        Set user for starting the servers. Default `whoami`.],
    [MEME_USER=$withval],
    [MEME_USER="$USER"])

# set contact for website
AC_ARG_WITH(
    contact,
    [  --with-contact=EMAIL    Set email of website administrator. Default `whoami`@`hostname`.],
    [CONTACT=$withval],
    [CONTACT="$MEME_USER@$HOSTNAME"])

# set website location
AC_ARG_WITH(
    webdir,
    [  --with-webdir=DIR       Set path to website part. Default $prefix/web.],
    [WEBDIR=$withval],
    [WEBDIR="\${prefix}/web"])

# set database location
AC_ARG_WITH(
    db,
    [  --with-db=DBDIR         Set path to databases. Default $prefix/db.],
    [MEME_DB=$withval],
    [MEME_DB="\${prefix}/db"])

# set logs location
AC_ARG_WITH(
    logs,
    [  --with-logs=LOGDIR      Set path to logs. Default $prefix/LOGS.],
    [MEME_LOGS=$withval],
    [MEME_LOGS="\${prefix}/LOGS"])

# set meme port
AC_ARG_WITH(
    meme-port,
    [  --with-meme-port=PORT   Set port to use for meme. Default 3012.],
    [MEME_PORT=$withval],
    [MEME_PORT="3012"])

# set mast port
AC_ARG_WITH(
    mast-port,
    [  --with-mast-port=PORT   Set port to use for meme. Default 3013.],
    [MAST_PORT=$withval],
    [MAST_PORT="3013"])

# set number of processors
AC_ARG_WITH(
    procs,
    [  --with-procs=PROCS      Set number of processors to use for parallel jobs. Default 2.],
    [NPROCS=$withval],
    [NPROCS="2"])

# set wall time limit 
AC_ARG_WITH(
    maxtime,
    [  --with-maxtime=TIME     Set wall time limit (sec) a job may use. Default 7200.],
    [MAXTIME=$withval],
    [MAXTIME="7200"])

# set queue name 
AC_ARG_WITH(
    qname,
    [  --with-qname=NAME       Set queue name (needed only for PBS scheduler). Default is empty string.],
    [QNAME=$withval],
    [QNAME=""])


CFLAGS="$CFLAGS $DEBUG_FLAGS $OPT_FLAGS"
MPIFLAGS="-DPARALLEL $ac_mpicc_flags"
MEME_DIRECTORY="\${prefix}"

AC_SUBST(MEME_DIRECTORY)
AC_SUBST(MACHID)
AC_SUBST(MEME_PORT)
AC_SUBST(MAST_PORT)
AC_SUBST(DEBUG_FLAGS)
AC_SUBST(OPT_FLAGS)
AC_SUBST(LIBS)
AC_SUBST(MPICC)
AC_SUBST(MPIINC)
AC_SUBST(MPI_CMD)
AC_SUBST(MPILIBDIR)
AC_SUBST(MPIFLAGS)
AC_SUBST(MPITYPE)
AC_SUBST(VERSION)
AC_SUBST(LD)
AC_SUBST(LAMBIN)
AC_SUBST(URL)
AC_SUBST(CONTACT)
AC_SUBST(MEME_USER)
AC_SUBST(MAIL)
AC_SUBST(SENDMAIL)
AC_SUBST(HOSTNAME)
AC_SUBST(SERVER)
AC_SUBST(WEBDIR)
AC_SUBST(MEME_DB)
AC_SUBST(MEME_LOGS)
AC_SUBST(QTYPE)
AC_SUBST(NPROCS)
AC_SUBST(MAXTIME)
AC_SUBST(QNAME)

AC_CONFIG_FILES([Makefile etc/Makefile scripts/Makefile src/Makefile src/parallel/Makefile tests/Makefile website/Makefile website/cgi-bin/Makefile website/html/Makefile website/scripts/Makefile website/examples/Makefile website/images/Makefile])
AC_OUTPUT

echo \
"
        ========================
        Configuration parameters
        ========================

  Install path:               ${prefix}
  Install UID:                ${MEME_USER}
  Version:                    ${VERSION}
  C compiler:                 ${CC}
  C compiler flags:           ${CFLAGS} 
  Linker:                     ${LD}
  Special Libs:               ${LIBS}
  MPICC                       ${MPICC}
  MPIINC                      ${MPIINC}
  MPILIBDIR                   ${MPILIBDIR}
  MPIFLAGS                    ${MPIFLAGS}
  MPIRUN                      ${MPI_CMD}
  URL:                        ${URL}
  PERL:                       ${PERL}
  MEME_LOGS:                  ${MEME_LOGS}

  Run the following commands to compile, test and install meme:
        make   
        make test
        make install
"
