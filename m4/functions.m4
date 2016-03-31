# $Id: functions.m4,v 1.1 2005/10/06 20:06:01 nadya Exp $
#
# $Log: functions.m4,v $
# Revision 1.1  2005/10/06 20:06:01  nadya
# add definitions for configure macros that are missing from autoconf v < 2.54
# will need to run bootstrap on install hostbefore running configure.
#
# This file is automatically added by bootstrap if needed

# This macro is only needed in autoconf <= 2.54.  Newer versions of autoconf
# have this macro built-in. This macro is copied from the file functions.m4
# from the autoconf-2.59 distribution.

# _AC_FUNC_REALLOC_IF(IF-WORKS, IF-NOT)
# ------------------------------------- 
# If `realloc (0, 0)' properly handled, run IF-WORKS, otherwise, IF-NOT.
AC_DEFUN([_AC_FUNC_REALLOC_IF],
[AC_REQUIRE([AC_HEADER_STDC])dnl
AC_CHECK_HEADERS(stdlib.h)
AC_CACHE_CHECK([for GNU libc compatible realloc], ac_cv_func_realloc_0_nonnull,
[AC_RUN_IFELSE(
[AC_LANG_PROGRAM(
[[#if STDC_HEADERS || HAVE_STDLIB_H
# include <stdlib.h>
#else
char *realloc ();
#endif
]],
                 [exit (realloc (0, 0) ? 0 : 1);])],
               [ac_cv_func_realloc_0_nonnull=yes],
               [ac_cv_func_realloc_0_nonnull=no],
               [ac_cv_func_realloc_0_nonnull=no])]) 
AS_IF([test $ac_cv_func_realloc_0_nonnull = yes], [$1], [$2])
])# AC_FUNC_REALLOC

# AC_FUNC_REALLOC
# ---------------
# Report whether `realloc (0, 0)' properly handled, and replace realloc if
# needed.
AN_FUNCTION([realloc], [AC_FUNC_REALLOC])
AC_DEFUN([AC_FUNC_REALLOC],
[_AC_FUNC_REALLOC_IF(
  [AC_DEFINE([HAVE_REALLOC], 1,
             [Define to 1 if your system has a GNU libc compatible `realloc'
              function, and to 0 otherwise.])],
  [AC_DEFINE([HAVE_REALLOC], 0)
   AC_LIBOBJ([realloc])
   AC_DEFINE([realloc], [rpl_realloc],
      [Define to rpl_realloc if the replacement function should be used.])])
])# AC_FUNC_REALLOC


