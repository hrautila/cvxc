
AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(srcdir)/include -I$(top_srcdir)/src/include
AM_CFLAGS = -O3 -march=native -ffast-math -fomit-frame-pointer -Wall

LIBFILES = \
	lib/conelp.c \
	lib/cp.c \
	lib/cpl.c \
	lib/cvxm.c \
	lib/cvxmio.c \
	lib/dims.c \
	lib/gp.c \
	lib/gpindex.c \
	lib/index.c \
	lib/json.c \
	lib/kkt.c \
	lib/kktldl.c \
	lib/kktldl2.c \
	lib/mgrp.c \
	lib/norm.c \
	lib/pack.c \
	lib/print.c \
	lib/scale.c \
	lib/scaling.c \
	lib/sfunc.c \
	lib/step.c

lib_LTLIBRARIES = libcvxc.la

libcvxc_la_SOURCES = $(LIBFILES)
libcvxc_la_LDFLAGS = -version-info @SOVERSION@

include_HEADERS = include/cvxc.h

noinst_HEADERS = \
	include/internal.h \
	include/cvxm.h

LIBCVXC = libcvxc.la

bin_PROGRAMS = cvxsolver
cvxsolver_SOURCES = bin/solvermain.c bin/input.c bin/shared.c bin/solve.c
cvxsolver_DEPENDENCIES = libcvxc.la
cvxsolver_CFLAGS = -I$(srcdir)/bin
cvxsolver_LDADD =  -L$(top_builddir)/src/.libs -lcvxc -larmasd -ldl
cvxsolver_LDFLAGS = -static




