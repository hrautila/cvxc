
INCLUDE = -I$(top_srcdir)/src/include -I$(srcdir)
AM_CFLAGS = -O3 -march=native -Wall $(INCLUDE)

LIBCVXC = $(top_builddir)/src/libcvxc.la

check_PROGRAMS = test_lp test_conelp test_sdp test_socp test_cpl test_cp test_cp2

test_lp_SOURCES = test_lp.c helper.c
test_lp_LDADD   = $(LIBCVXC) -larmasd
test_conelp_SOURCES = test_conelp.c helper.c
test_conelp_LDADD   = $(LIBCVXC) -larmasd
test_sdp_SOURCES = test_sdp.c helper.c
test_sdp_LDADD   = $(LIBCVXC) -larmasd
test_socp_SOURCES = test_socp.c helper.c
test_socp_LDADD   = $(LIBCVXC) -larmasd

test_cpl_SOURCES = test_cpl.c helper.c
test_cpl_LDADD   = $(LIBCVXC) -larmasd
test_cp_SOURCES = test_cp.c helper.c
test_cp_LDADD   = $(LIBCVXC) -larmasd

test_cp2_SOURCES = test_cp2.c helper.c
test_cp2_LDADD   = $(LIBCVXC) -larmasd


TESTS = check_PROGRAMS
