
##NETLIB = ../armas/src/.libs/libarmasd.so
NETLIB = ../armas/src/.libs/libarmasd.a

CFLAGS = -g -O0 -I../armas/src -Wall -pthread

CVXOBJ = dims.o jfunc.o sfunc.o scale.o conelp.o kkt.o misc.o
MATOBJ = cvxm.o cvxmio.o

LIBOBJ = $(CVXOBJ) $(MATOBJ)
OBJS   = test_lp.o test_conelp.o test_socp.o test_sdp.o

$(LIBOBJ): convex.h cvxm.h
$(OBJS): convex.h cvxm.h

libcvx.a : $(LIBOBJ)
	$(AR) rs $@ $(LIBOBJ)

test_lp: test_lp.o libcvx.a
	$(CC) $(CFLAGS) -o $@ test_lp.o libcvx.a $(NETLIB) -lm
test_socp: test_socp.o libcvx.a
	$(CC) $(CFLAGS) -o $@ test_socp.o libcvx.a $(NETLIB) -lm
test_sdp: test_sdp.o libcvx.a
	$(CC) $(CFLAGS) -o $@ test_sdp.o libcvx.a $(NETLIB) -lm
test_conelp: test_conelp.o libcvx.a
	$(CC) $(CFLAGS) -o $@ test_conelp.o libcvx.a $(NETLIB) -lm

test_init: test_init.o libcvx.a
	$(CC) $(CFLAGS) -o $@ test_init.o libcvx.a $(NETLIB) -lm

test_dims: test_dims.c cvxm.c dims.c convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_dims.c cvxm.c dims.c $(NETLIB) -lm

test_jfunc: test_jfunc.c jfunc.c cvxm.c convex.h cvxm.h dims.c cvxmio.c
	$(CC) $(CFLAGS) -o $@ test_jfunc.c jfunc.c cvxm.c dims.c cvxmio.c $(NETLIB) -lm

test_sfunc: test_sfunc.c sfunc.c cvxm.c dims.c convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_sfunc.c sfunc.c cvxm.c dims.c $(NETLIB) -lm

test_scale: test_scale.o libcvx.a convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_scale.o libcvx.a $(NETLIB) -lm


faa: faa.c
	$(CC) $(CFLAGS) -o $@ faa.c $(NETLIB) -lm
