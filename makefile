
##NETLIB = ../armas/src/.libs/libarmasd.so
NETLIB = ../armas/src/.libs/libarmasd.a

CFLAGS = -g -O0 -I../armas/src -Wall -pthread

test_dims: test_dims.c cvxm.c dims.c convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_dims.c cvxm.c dims.c $(NETLIB) -lm

test_jfunc: test_jfunc.c jfunc.c cvxm.c convex.h cvxm.h dims.c cvxmio.c
	$(CC) $(CFLAGS) -o $@ test_jfunc.c jfunc.c cvxm.c dims.c cvxmio.c $(NETLIB) -lm

test_sfunc: test_sfunc.c sfunc.c cvxm.c dims.c convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_sfunc.c sfunc.c cvxm.c dims.c $(NETLIB) -lm

test_scale: test_scale.o scale.o sfunc.o jfunc.o cvxm.o dims.o cvxmio.o convex.h cvxm.h
	$(CC) $(CFLAGS) -o $@ test_scale.o scale.o sfunc.o jfunc.o cvxm.o dims.o cvxmio.o $(NETLIB) -lm


