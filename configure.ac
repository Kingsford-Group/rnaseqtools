#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT(rnaseqtools, 1.0.2, mingfu.shao@cs.cmu.edu)
AM_INIT_AUTOMAKE()
AC_CONFIG_SRCDIR([gtfalign/config.h])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET
AC_PROG_RANLIB

# Checks for libraries.

# Checks for header files.
AC_CHECK_HEADERS([stdint.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_CHECK_FUNCS([sqrt])

AC_CONFIG_FILES([Makefile
				 gtfalign/Makefile
                 gtfcuff/Makefile
                 gtfformat/Makefile
                 gtfmerge/Makefile
                 gtfquant/Makefile
                 gtfsimulator/Makefile
                 lib/Makefile
                 lib/gtf/Makefile
                 lib/util/Makefile])
AC_OUTPUT
