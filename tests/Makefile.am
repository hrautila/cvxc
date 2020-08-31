
INCLUDE = -I$(top_srcdir)/src/include -I$(srcdir)
AM_CFLAGS = -O3 -march=native -Wall $(INCLUDE)

LIBCVXC = $(top_builddir)/src/libcvxc.la

check_PROGRAMS = test_lp test_conelp

test_lp_SOURCES = test_lp.c
test_lp_LDADD   = $(LIBCVXC) -larmasd
test_conelp_SOURCES = test_conelp.c
test_conelp_LDADD   = $(LIBCVXC) -larmasd

TESTS = check_PROGRAMS
