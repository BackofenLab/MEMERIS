# $Id: headers.m4,v 1.1 2005/10/06 20:06:01 nadya Exp $
#
# $Log: headers.m4,v $
# Revision 1.1  2005/10/06 20:06:01  nadya
# add definitions for configure macros that are missing from autoconf v < 2.54
# will need to run bootstrap on install hostbefore running configure.
#
# This file is automatically added by bootstrap if needed

# This macro is only needed in autoconf <= 2.54.  Newer versions of autoconf
# have this macro built-in. This macro is copied from the file headers.m4
# from the autoconf-2.59 distribution.

AC_DEFUN([AC_HEADER_STDBOOL],
  [AC_CACHE_CHECK([for stdbool.h that conforms to C99],
     [ac_cv_header_stdbool_h],
     [AC_TRY_COMPILE(
        [
          #include <stdbool.h>
          #ifndef bool
           "error: bool is not defined"
          #endif
          #ifndef false
           "error: false is not defined"
          #endif
          #if false
           "error: false is not 0"
          #endif
          #ifndef true
           "error: false is not defined"
          #endif
          #if true != 1
           "error: true is not 1"
          #endif
          #ifndef __bool_true_false_are_defined
           "error: __bool_true_false_are_defined is not defined"
          #endif

          struct s { _Bool s: 1; _Bool t; } s;

          char a[true == 1 ? 1 : -1];
          char b[false == 0 ? 1 : -1];
          char c[__bool_true_false_are_defined == 1 ? 1 : -1];
          char d[(bool) -0.5 == true ? 1 : -1];
          bool e = &s;
          char f[(_Bool) -0.0 == false ? 1 : -1];
          char g[true];
          char h[sizeof (_Bool)];
          char i[sizeof s.t];
        ],
        [ return !a + !b + !c + !d + !e + !f + !g + !h + !i; ],
        [ac_cv_header_stdbool_h=yes],
        [ac_cv_header_stdbool_h=no])])
   AC_CHECK_TYPES([_Bool])
   if test $ac_cv_header_stdbool_h = yes; then
     AC_DEFINE(HAVE_STDBOOL_H, 1, [Define to 1 if stdbool.h conforms to C99.])
   fi])
