dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.
dnl
AC_PREREQ([2.69])
AC_INIT([CvxOpt in C], [m4_esyscmd(./package-version.sh)],
   [harri.rautila@gmail.com], [libcvx], [])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/Makefile.am])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT

dnl shared library version number
AC_SUBST(SOVERSION, [m4_esyscmd(./abi-version.sh -libtool)])

dnl Checks for programs.
AC_PROG_CC

dnl Checks for libraries.
AC_CHECK_LIB([m], [fabs])
AC_CHECK_LIB([armasd], [armas_version], [libarmas_found=yes], [])
AS_IF([test "x$libarmas_found" != "xyes"],
      [AC_MSG_ERROR([Unable to find libarmasd files. Define LDFLAGS=-L/path/includes on invocation.])])

dnl Checks for header files.
AC_CHECK_HEADER(armas/ddense.h, [armas_dense_found=yes])
AS_IF([test "x$armas_dense_found" != "xyes"],
      [AC_MSG_ERROR([Unable to find libarmas header files. Define CPPFLAGS=-I/path/includes on invocation.])])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT64_T
AC_C_RESTRICT
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT64_T

dnl Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday sqrt])

AC_CONFIG_FILES([
   Makefile
   src/Makefile
])

AC_OUTPUT
