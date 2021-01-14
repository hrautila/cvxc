dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.
dnl
AC_PREREQ([2.69])
AC_INIT([CVXopt in C], [m4_esyscmd(./package-version.sh)],
   [harri.rautila@gmail.com], [cvxc], [https://github.com/hrautila/cvxc])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/Makefile.am])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])
LT_INIT

dnl shared library version number
AC_SUBST(SOVERSION, [m4_esyscmd(./abi-version.sh -libtool)])

dnl Checks for programs.
AC_PROG_CC

AC_ARG_ENABLE(float32,
   AC_HELP_STRING([--enable-float32],[enable single precission float (default is no)]),
      [enable_float32=$enableval], [enable_float32=no])

AS_IF([test "$enable_float32" = "no"], [libarmas=armasd], [libarmas=armass])
AS_IF([test "$enable_float32" = "no"], [armashdr="ddense.h"], [armashdr="sdense.h"])

dnl Checks for libraries.
AC_CHECK_LIB([m], [fabs])
AC_CHECK_LIB([$libarmas], [armas_version], [libarmas_found=yes], [])
AS_IF([test "x$libarmas_found" != "xyes"],
      [AC_MSG_ERROR([Unable to find libarmas files. Define LDFLAGS=-L/path/includes on invocation.])])

dnl Checks for header files.
AC_CHECK_HEADER("armas/$armashdr", [armas_dense_found=yes])
AS_IF([test "x$armas_dense_found" != "xyes"],
      [AC_MSG_ERROR([Unable to find libarmas header files. Define CPPFLAGS=-I/path/includes on invocation.])])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT64_T
AC_CHECK_HEADER_STDBOOL
AC_TYPE_UINT64_T

dnl Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday sqrt])

AM_CONDITIONAL([ENABLE_FLOAT32],    [test "$enable_float32" = "yes"])

AC_CONFIG_FILES([
   Makefile
   src/Makefile
   tests/Makefile
])

AC_OUTPUT
